/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.cobi.genetsimulator.entity;

/**
 *
 * @author mxli
 */
public class PValueSet {
     int count;
     double proportion;

    public PValueSet(int count, double proportion) {
        this.count = count;
        this.proportion = proportion;
    }

    public int getCount() {
        return count;
    }

    public void setCount(int count) {
        this.count = count;
    }

    public double getProportion() {
        return proportion;
    }

    public void setProportion(double proportion) {
        this.proportion = proportion;
    }
     
}
