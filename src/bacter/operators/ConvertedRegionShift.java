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

package bacter.operators;

import bacter.Conversion;
import beast.core.Description;
import beast.core.Input;
import beast.util.Randomizer;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Operator which moves the alignment region affected "
        + "by a randomly-selected conversion event.")
public class ConvertedRegionShift extends ACGOperator {

    public Input<Double> apertureSizeInput = new Input<>(
            "apertureSize",
            "Relative size (with respect to alignment size) of aperture "
                    + "within which new location of region edge is chosen "
                    + "uniformly. (Default 0.01, ie. 1%)", 0.01);

    public ConvertedRegionShift() { }

    @Override
    public double proposal() {
        
        if (acg.getTotalConvCount()<1 || acg.wholeLocusModeOn())
            return Double.NEGATIVE_INFINITY;

        Conversion conv = chooseConversion();
        
        int radius = (int)Math.round(conv.getLocus().getSiteCount()
            *apertureSizeInput.get())/2;

        int delta = Randomizer.nextInt(radius*2 + 1) - radius;

        int newStart = conv.getStartSite()+delta;
        int newEnd = conv.getEndSite()+delta;

        if (!acg.circularGenomeModeOn()) {                                                      //todo: check adjustment (circular genome)
            if (newEnd > conv.getLocus().getSiteCount() - 1)
                return Double.NEGATIVE_INFINITY;
            if (newStart < 0)
                return Double.NEGATIVE_INFINITY;
        } else {
            if (newStart < 0 || newStart > acg.getTotalConvertibleSequenceLength()-1) {
                newStart = newStart < 0 ? acg.getTotalConvertibleSequenceLength() + newStart : newStart-acg.getTotalConvertibleSequenceLength();
            }
            if (newEnd < 0 || newEnd > acg.getTotalConvertibleSequenceLength()-1) {
                newEnd = newEnd < 0 ? acg.getTotalConvertibleSequenceLength() + newEnd : newEnd-acg.getTotalConvertibleSequenceLength() ;
            }
        }

        conv.setStartSite(newStart);
        conv.setEndSite(newEnd);

        assert !acg.isInvalid() : "CRS produced invalid state.";
        
        return 0.0;
    }
    
}
