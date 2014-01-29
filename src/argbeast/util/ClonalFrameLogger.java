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

package argbeast.util;

import argbeast.RecombinationGraph;
import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Loggable;
import java.io.PrintStream;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Logs clonal frame corresponding to recombination graph.")
public class ClonalFrameLogger extends CalculationNode implements Loggable {

    public Input<RecombinationGraph> argInput = new Input<RecombinationGraph>(
            "arg", "Recombination graph whose clonal frame you want to log.",
            Validate.REQUIRED);

    @Override
    public void initAndValidate() { }
    
    @Override
    public void init(PrintStream out) throws Exception {
       argInput.get().init(out);
    }

    @Override
    public void log(int nSample, PrintStream out) {
        RecombinationGraph arg = argInput.get();

        out.print("tree STATE_" + nSample + " = ");
        out.print(arg.getRoot().toSortedNewick(new int[1], false) + ";");
    }

    @Override
    public void close(PrintStream out) {
        argInput.get().close(out);
    }
    
}
