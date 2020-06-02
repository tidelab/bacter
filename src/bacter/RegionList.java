/*
 * Copyright (C) 2014 Tim Vaughan <tgvaughan@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package bacter;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

import java.util.*;
import java.util.function.Consumer;

/**
 * This class is used to maintain a list of marginal tree regions
 * corresponding to a given ACG.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class RegionList {

    private final List<Region> regions;
    private Locus locus;
    private boolean dirty;

    /**
     * Ancestral conversion graph this list belongs to.
     */
    private final ConversionGraph acg;

    /**
     * Construct a new region list for the given ACG.  There should only
     * be one of these objects per ACG object, created during the ACG
     * initAndValidate().
     *
     * @param acg Conversion graph from which to compute region list.
     * @param locus Locus with which this region list is associated.
     */
    public RegionList(ConversionGraph acg, Locus locus) {
        this.acg = acg;
        this.locus = locus;
        regions = new ArrayList<>();
        dirty = true;
    }

    /**
     * Obtain list of contiguous regions having fixed marginal trees.
     * 
     * @return region list
     */
    public List<Region> getRegions() {
        updateRegionList();

        return regions;
    }

    /**
     * Retrieve number of regions in region list.
     * 
     * @return region count
     */
    public int getRegionCount() {
        updateRegionList();

        return regions.size();
    }

    /**
     * Mark the region list as dirty.
     */
    public void makeDirty() {
        dirty = true;
    }
   
    /**
     * Assemble list of regions of contiguous sites that possess a single
     * marginal tree.
     */
    public void updateRegionList() {
        if (!dirty)
            return;

        regions.clear();

        AffectedSiteList affectedSiteList = new AffectedSiteList(acg);

        /* Assemble lists of conversions ordered by start and end sites.
        Note that these are COPIES of the conversion objects attached
        to the ACG. This ensures that subsequent modifications of these
        objects won't break our contract with the HashSet<Conversion>
        objects in the likelihood code.
        */
        List<Conversion> convOrderedByStart = new ArrayList<>();

        acg.getConversions(locus).forEach(conversion -> {
            if (affectedSiteList.affectedSiteCount.get(conversion)>0)
                convOrderedByStart.add(conversion.getCopy());
        });
        convOrderedByStart.sort(Comparator.comparingInt((Conversion o) -> o.startSite));

        List<Conversion> convOrderedByEnd = new ArrayList<>();
        convOrderedByEnd.addAll(convOrderedByStart);
        convOrderedByEnd.sort(Comparator.comparingInt((Conversion o) -> o.endSite));

        Set<Conversion> activeConversions = Sets.newHashSet();

        //todo: revise and check adjustment (circular genome)
        int lastBoundary = 0;
        boolean convOverlap = false, firstStep = false, noConv = true;
        if (acg.circularGenomeModeOn() && !convOrderedByStart.isEmpty()) {
            lastBoundary = Math.max(convOrderedByEnd.get(convOrderedByEnd.size() - 1).getEndSite() + 1, convOrderedByStart.get(convOrderedByStart.size() - 1).getStartSite());
            convOverlap = (convOrderedByEnd.get(0).getEndSite() + 1) < (convOrderedByStart.get(0).getStartSite());
            noConv = false;
            firstStep = !convOverlap;
        }

        while (!convOrderedByStart.isEmpty() || !convOrderedByEnd.isEmpty()) {

            int nextStart;
            if (!convOrderedByStart.isEmpty())
                nextStart = convOrderedByStart.get(0).getStartSite();
            else
                nextStart = Integer.MAX_VALUE;

            int nextEnd;
            if (!convOrderedByEnd.isEmpty())
                nextEnd = convOrderedByEnd.get(0).getEndSite() + 1;
            else
                nextEnd = Integer.MAX_VALUE;

            int nextBoundary = Math.min(nextStart, nextEnd);
            if (convOverlap) {
                activeConversions.add(convOrderedByEnd.get(0));
                Region region = new Region(lastBoundary, nextBoundary, activeConversions, acg.getTotalConvertibleSequenceLength());
                regions.add(region);
                activeConversions.remove(convOrderedByEnd.get(0));
                convOrderedByEnd.remove(0);
            }
            if (nextBoundary > lastBoundary || firstStep) {
                if (nextBoundary == 0) {
                    nextBoundary = locus.getSiteCount();
                }
                Region region = new Region(lastBoundary, nextBoundary, activeConversions, acg.getTotalConvertibleSequenceLength());
                regions.add(region);
                firstStep = false;
            }

            if (nextStart < nextEnd) {
                activeConversions.add(convOrderedByStart.get(0));
                convOrderedByStart.remove(0);
                lastBoundary = nextStart;
            } else if (!convOverlap) {
                activeConversions.remove(convOrderedByEnd.get(0));
                convOrderedByEnd.remove(0);
                lastBoundary = nextEnd;
            } else {
                lastBoundary = nextBoundary;
                convOverlap = false;
                firstStep = true;
            }
        }

        if ((lastBoundary < locus.getSiteCount()) && noConv) {
            Region region = new Region(lastBoundary, locus.getSiteCount(), new HashSet<>(), acg.getTotalConvertibleSequenceLength());
            //Region region = new Region(lastBoundary, locus.getSiteCount(), new HashSet<>());
            regions.add(region);
        }

        dirty = false;
    }
}
