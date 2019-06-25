package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import com.google.common.annotations.VisibleForTesting;
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

    private int weightThresholdToUse = DEFAULT_OUTGOING_JT_EVIDENCE_THRESHOLD_TO_BELEIVE;

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

    @VisibleForTesting
    public void setWeightThresholdToUse(final int outgoingWeight) {
        weightThresholdToUse = outgoingWeight;
    }

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
                    outgoingEdges = graph.outgoingEdgesOf(vertexToExtend);
                }
                // Cache the chain result
                contiguousSequences.put(pathToExtend.getLastVertex(), chain);

            } else {
                // We have already expanded this part of the graph, use it.
                vertexToExtend = graph.getEdgeTarget(chain.get(chain.size() - 1));
            }
            // vertexToExtend is necessarily at the next "interesting" point

            Collection<E> outgoingEdges = graph.outgoingEdgesOf(vertexToExtend);

            // In the event we have a junction tree on top of a vertex with outDegree > 1, we add this first before we traverse paths
            if ( experimentalReadThreadingGraph.getJunctionTreeForNode((MultiDeBruijnVertex) vertexToExtend) != null) { //TODO make the condition for this actually based on the relevant junction tree
                // TODO chain can be null but we still need to inherit a thing, probably happens whenever we pick up a tree.
                pathToExtend.addJunctionTree(experimentalReadThreadingGraph.getJunctionTreeForNode((MultiDeBruijnVertex) vertexToExtend));
            }

            // If we are at a reference end then we close out the path TODO this isn't adequate for non-unique reference sinks
            if (sinks.contains(vertexToExtend)) {
                //TODO this will probably be resolved using a junction tree on that node and treating it as an edge to extend
                //todo the proposal here would be to check if there is an active tree left for us at this point and if so keep going
                if (chain.isEmpty()) {
                    result.add(pathToExtend);
                } else {
                    result.add(new RTBesthaplotype<>(pathToExtend, chain, 0));
                }

            // We must be at a point where the path diverges, use junction trees to resolve if possible
            } else {
                if (outgoingEdges.size() > 1) {
                    List<RTBesthaplotype<V, E>> jTPaths = pathToExtend.getApplicableNextEdgesBasedOnJunctionTrees(chain, weightThresholdToUse);
                    if (jTPaths.isEmpty()) {
                        // Standard behavior from the old KBestHaplotypeFinder
                        int totalOutgoingMultiplicity = 0;
                        for (final BaseEdge edge : outgoingEdges) {
                            totalOutgoingMultiplicity += edge.getMultiplicity();
                        }

                        // Add all valid edges to the graph
                        for (final E edge : outgoingEdges) {
                            // Don't traverse an edge if it only has reference evidence supporting it (unless there is no other evidence whatsoever)
                            if (totalOutgoingMultiplicity != 0 && edge.getMultiplicity() != 0) {
                                List<E> chainCopy = new ArrayList<>(chain);
                                chainCopy.add(edge);
                                queue.add(new RTBesthaplotype<>(pathToExtend, chainCopy, edge.getMultiplicity(), totalOutgoingMultiplicity));
                            }
                        }
                    } else {
                        queue.addAll(jTPaths);
                    }

                // Otherwise there is only one node forwards so we take that
                } else {
                    List<E> chainCopy = new ArrayList<>(chain);
                    chainCopy.add(outgoingEdges.iterator().next());
                    queue.add(new RTBesthaplotype<>(pathToExtend, chainCopy, 0));
                }
            }
        }

        return result.stream().map(n -> (KBestHaplotype<V, E>) n).collect(Collectors.toList());
    }

    /**
     * Removes edges that produces cycles and also dead vertices that do not lead to any sink vertex.
     * @return never {@code null}.
     */
    @Override
    protected BaseGraph<V, E> removeCyclesAndVerticesThatDontLeadToSinks(final BaseGraph<V, E> original, final Collection<V> sources, final Set<V> sinks) {
        return original;
    }

}
