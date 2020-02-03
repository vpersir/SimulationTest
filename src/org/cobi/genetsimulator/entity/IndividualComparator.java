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
public class IndividualComparator implements Comparator<Individual> {

    public int compare(Individual arg0, Individual arg1) {
        return Double.compare(arg0.getMainTrait()[0], arg1.getMainTrait()[0]);
    }
}
