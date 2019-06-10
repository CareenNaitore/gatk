package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.apache.commons.lang3.mutable.MutableInt;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ExperimentalReadThreadingGraph;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.MultiDeBruijnVertex;
import org.broadinstitute.hellbender.utils.Utils;
import org.jgrapht.alg.CycleDetector;

import java.util.*;
import java.util.stream.Collectors;

public class ExperimentalKBestHaplotypeFinder<V extends BaseVertex, E extends BaseEdge> extends  KBestHaplotypeFinder<V, E> {
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

//    /**
//     * Implement Dijkstra's algorithm as described in https://en.wikipedia.org/wiki/K_shortest_path_routing
//     */
//    @Override
//    public List<KBestHaplotype<V, E>> findBestHaplotypes(final int maxNumberOfHaplotypes) {
//        final List<KBestHaplotype<V, E>> result = new ArrayList<>();
//        asdfasdf
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
