/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.genetsimulator.entity;

import java.util.Comparator;

/**
 *
 * @author mxli
 */
public class PValueWeightWComparator1 implements Comparator <PValueWeight>{

    public int compare(PValueWeight arg0, PValueWeight arg1) { 
        return Double.compare(arg1.var, arg0.var);
    }
}
