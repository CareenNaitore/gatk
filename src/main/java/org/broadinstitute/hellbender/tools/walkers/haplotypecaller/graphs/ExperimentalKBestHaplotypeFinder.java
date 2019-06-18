package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.apache.commons.lang3.mutable.MutableInt;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ExperimentalReadThreadingGraph;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.MultiDeBruijnVertex;
import org.broadinstitute.hellbender.utils.Utils;
import org.jgrapht.alg.CycleDetector;

import java.util.*;
import java.util.stream.Collectors;

public class ExperimentalKBestHaplotypeFinder<V extends BaseVertex, E extends BaseEdge> extends  KBestHaplotypeFinder<V, E> {
    // TODO PARAMETERETERIZE ME
    public static final int DEFAULT_OUTGOING_JT_EVIDENCE_THRESHOLD_TO_BELEIVE = 3;

    ExperimentalReadThreadingGraph experimentalReadThreadingGraph;

    public ExperimentalKBestHaplotypeFinder(BaseGraph<V, E> graph, Set<V> sources, Set<V> sinks) {
        super(graph, sources, sinks);
        if (graph instanceof ExperimentalReadThreadingGraph) {
            experimentalReadThreadingGraph = (ExperimentalReadThreadingGraph) graph;
        } else {
            throw new RuntimeException("ExperimentalKBesthaplotypeFinder requires an ExperimentalReadThreadingGraph be provided");
        }
    }

    /**
     * Constructor for the special case of a single source and sink
     */
    public ExperimentalKBestHaplotypeFinder(final BaseGraph<V, E> graph, final V source, final V sink) {
        this(graph, Collections.singleton(source), Collections.singleton(sink));
    }

    /**
     * Constructor for the default case of all sources and sinks
     */
    public ExperimentalKBestHaplotypeFinder(final BaseGraph<V, E> graph) {
        this(graph, graph.getReferenceSourceVertex(), graph.getReferenceSinkVertex());
    }

    // We want to accept graphs with cycles at this stage if we think they are still resolvable
    @Override
    protected BaseGraph<V, E> removeCyclesIfNecessary(BaseGraph<V, E> graph, Set<V> sources, Set<V> sinks) {
        return graph;
    }

    Map<V, List<E>> contiguousSequences = new HashMap<>();

    @Override
    @SuppressWarnings({"unchecked"})
    public List<KBestHaplotype<V, E>> findBestHaplotypes(final int maxNumberOfHaplotypes) {
        final List<RTBesthaplotype<V, E>> result = new ArrayList<>();
        final PriorityQueue<RTBesthaplotype<V, E>> queue = new PriorityQueue<>(Comparator.comparingDouble(KBestHaplotype<V, E>::score).reversed());
        sources.forEach(source -> queue.add(new RTBesthaplotype<>(source, graph)));

        //TODO this optimization has been disabled for now to simplify the code, this will probably need to be reenabled later
//        final Map<V, MutableInt> vertexCounts = graph.vertexSet().stream()
//                .collect(Collectors.toMap(v -> v, v -> new MutableInt(0)));

        // Iterate over paths in the queue
        while (!queue.isEmpty() && result.size() < maxNumberOfHaplotypes) {
            final RTBesthaplotype<V, E> pathToExtend = queue.poll();
            V vertexToExtend = pathToExtend.getLastVertex();
            final List<E> chain = contiguousSequences.computeIfAbsent(vertexToExtend, k -> new ArrayList<>());

            // If we have a cached chain use its target vertex
            if (chain.isEmpty()){
                Set<E> outgoingEdges = graph.outgoingEdgesOf(vertexToExtend);
                // Keep going until we reach a fork, reference sink, or fork
                while (!sinks.contains(vertexToExtend) && outgoingEdges.size() == 1 &&
                        experimentalReadThreadingGraph.getJunctionTreeForNode((MultiDeBruijnVertex) vertexToExtend) == null) {
                    E edge = outgoingEdges.iterator().next();
                    chain.add(edge);
                    vertexToExtend = graph.getEdgeTarget(edge);
                }
                // Cache the chain result
                contiguousSequences.put(pathToExtend.getLastVertex(), chain);

            } else {
                // We have already expanded this part of the graph, use it.
                vertexToExtend = graph.getEdgeTarget(chain.get(chain.size() - 1));
            }
            // vertexToExtend is necessarily at the next "interesting" point

            Collection<E> outgoingEdges = graph.outgoingEdgesOf(vertexToExtend);

            // if we are at a reference end then we close out the path TODO this isn't adequate for non-unique reference sinks
            if (sinks.contains(vertexToExtend)) {
                //TODO this will probably be resolved using a junction tree on that node and treating it as an edge to extend
                //todo the proposal here would be to check if there is an active tree left for us at this point and if so keep going
                result.add(pathToExtend);

            // We are at the vertex before one with an inDegree > 1, then we extend the path and add the tree to the list
            } else if ( experimentalReadThreadingGraph.getJunctionTreeForNode((MultiDeBruijnVertex) vertexToExtend) != null) { //TODO make the condition for this actually based on the relevant junction tree
                queue.add(new RTBesthaplotype<>(pathToExtend, chain, 0, experimentalReadThreadingGraph.getJunctionTreeForNode((MultiDeBruijnVertex) vertexToExtend)));

            // We must be at a point where the path diverges, use junction trees to resolve if possible
            } else {
                List<RTBesthaplotype<V, E>> jTPaths = pathToExtend.getApplicableNextEdgesBasedOnJunctionTrees(outgoingEdges, DEFAULT_OUTGOING_JT_EVIDENCE_THRESHOLD_TO_BELEIVE);
                if (jTPaths.isEmpty()) {
                    // Standard behavior from the old KBestHaplotypeFinder
                    int totalOutgoingMultiplicity = 0;
                    for (final BaseEdge edge : outgoingEdges) {
                        totalOutgoingMultiplicity += edge.getMultiplicity();
                    }

                    for (final E edge : outgoingEdges) {
                        queue.add(new RTBesthaplotype<>(pathToExtend, edge, edge.getMultiplicity(), totalOutgoingMultiplicity));
                    }
                } else {
                    queue.addAll(jTPaths);
                }
            }
        }

        return result.stream().map(n -> (KBestHaplotype<V, E>) n).collect(Collectors.toList());
    }

//    /**
//     * Implement Dijkstra's algorithm as described in https://en.wikipedia.org/wiki/K_shortest_path_routing
//     */
//    @Override
//    public List<KBestHaplotype<V, E>> findBestHaplotypes(final int maxNumberOfHaplotypes) {
//        final List<KBestHaplotype<V, E>> result = new ArrayList<>();
//        final PriorityQueue<KBestHaplotype<V, E>> queue = new PriorityQueue<>(Comparator.comparingDouble(KBestHaplotype<V, E>::score).reversed());
//        sources.forEach(source -> queue.add(new KBestHaplotype<>(source, graph)));
//
//        final Map<V, MutableInt> vertexCounts = graph.vertexSet().stream()
//                .collect(Collectors.toMap(v -> v, v -> new MutableInt(0)));
//
//        while (!queue.isEmpty() && result.size() < maxNumberOfHaplotypes) {
//            final KBestHaplotype<V, E> pathToExtend = queue.poll();
//            final V vertexToExtend = pathToExtend.getLastVertex();
//            if (sinks.contains(vertexToExtend)) {
//                result.add(pathToExtend);
//            } else {
//                asdf
//                if (vertexCounts.get(vertexToExtend).getAndIncrement() < maxNumberOfHaplotypes) {
//                    final Set<E> outgoingEdges = graph.outgoingEdgesOf(vertexToExtend);
//                    int totalOutgoingMultiplicity = 0;
//                    for (final BaseEdge edge : outgoingEdges) {
//                        totalOutgoingMultiplicity += edge.getMultiplicity();
//                    }
//
//                    for (final E edge : outgoingEdges) {
//                        final V targetVertex = graph.getEdgeTarget(edge);
//                        queue.add(new KBestHaplotype<>(pathToExtend, edge, totalOutgoingMultiplicity));
//                    }
//                }
//            }
//        }
//        return result;
//    }

    /**
     * Removes edges that produces cycles and also dead vertices that do not lead to any sink vertex.
     * @return never {@code null}.
     */
    @Override
    protected BaseGraph<V, E> removeCyclesAndVerticesThatDontLeadToSinks(final BaseGraph<V, E> original, final Collection<V> sources, final Set<V> sinks) {
        return original;
    }

}
