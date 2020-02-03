/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.genetsimulator.controller;

import java.io.Serializable;

/**
 *
 * @author mxli
 */
public class LDPair implements Serializable {

    public int i1;
    public int i2;
    public float ld;

    public LDPair(int index1, int index2, float ld) {
        this.i1 = index1;
        this.i2 = index2;
        this.ld = ld;
    }
}
