package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ExperimentalReadThreadingGraph;

import java.util.*;
import java.util.stream.Collectors;

public class RTBesthaplotype<T extends BaseVertex, E extends BaseEdge> extends KBestHaplotype<T, E> {
    // TODO these nodes probably shouldn't be exposed and held externally methinks
    private List<ExperimentalReadThreadingGraph.ThreadingNode> activeNodes;

    public RTBesthaplotype(final RTBesthaplotype p, final List<E> edgesToExtend, final double edgePenalty, final ExperimentalReadThreadingGraph.ThreadingTree treeToAdd) {
        super(p, edgesToExtend, edgePenalty);
        //this.threadingTreesToConsult = p.threadingTreesToConsult.cop;
        activeNodes = new ArrayList<ExperimentalReadThreadingGraph.ThreadingNode>(p.activeNodes);
        if (treeToAdd != null) {
            //threadingTreesToConsult.add(treeToAdd);
            activeNodes.add(treeToAdd.getRootNode());
        }
    }

    // Constructor to be used for internal calls from {@link #getApplicableNextEdgesBasedOnJunctionTrees()}
    public RTBesthaplotype(final RTBesthaplotype p, final List<E> chain, final int edgeMultiplicity, final int totalOutgoingMultiplicity) {
        super(p, chain, computeLogPenaltyScore( edgeMultiplicity, totalOutgoingMultiplicity));
        activeNodes = new ArrayList<ExperimentalReadThreadingGraph.ThreadingNode>(p.activeNodes);
        // Ensure that the relevant edge has been traversed
        takeEdge(chain.get(chain.size() - 1));
    }

    public RTBesthaplotype(final T initialVertex, final BaseGraph<T,E> graph) {
        super(initialVertex, graph);
        activeNodes = new ArrayList<>();
    }

    /**
     * This method is the primary engine for parsing
     */
    @SuppressWarnings({"unchecked"})
    public List<RTBesthaplotype<T, E>> getApplicableNextEdgesBasedOnJunctionTrees(final Collection<E> edgesAtVertex, final int weightThreshold) {
        List<RTBesthaplotype<T, E>> output = new ArrayList<>();
        ExperimentalReadThreadingGraph.ThreadingNode eldestTree = activeNodes.isEmpty() ? null : activeNodes.get(0);
        while (eldestTree != null) {
            //TODO this can be better
            int totalOut = 0;
            for ( ExperimentalReadThreadingGraph.ThreadingNode node : eldestTree.getChildrenNodes().values()) {
                totalOut += node.getCount();
            }
            // This right here is what handles dealing with weight thresholds
            if (totalOut >= weightThreshold) {
                //TODO add SOME sanity check to ensure that the vertex we stand on and the edges we are polling line up
                for (Map.Entry<MultiSampleEdge, ExperimentalReadThreadingGraph.ThreadingNode> childNode : eldestTree.getChildrenNodes().entrySet()) {
                    ExperimentalReadThreadingGraph.ThreadingNode child = childNode.getValue();
                    output.add(new RTBesthaplotype<>(this, Collections.singletonList((E) childNode.getKey()), child.getCount(), totalOut));
                }
                return output;

            // If there aren't enough outgoing nodes we just remove
            } else {
                activeNodes.remove(0);
                eldestTree = activeNodes.isEmpty() ? null : activeNodes.get(0);
            }
        }
        return output;
    }

    // method to handle incrementing all of the nodes in the tree
    private void takeEdge(final E edgeTaken) {
        activeNodes = activeNodes.stream().map(node -> {
            if (!node.getChildrenNodes().containsKey(edgeTaken)) {
                return null;
            }
            return node.getChildrenNodes().get(edgeTaken);
        }).filter(Objects::nonNull).collect(Collectors.toList());
    }

    public void addJunctionTree(final ExperimentalReadThreadingGraph.ThreadingTree junctionTreeForNode) {
        activeNodes.add(junctionTreeForNode.getRootNode());
    }
}
