/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.cobi.genetsimulator.controller;

import java.io.Serializable;
import java.util.Comparator;

/**
 *
 * @author mxli
 */
public class SparseLDMatrixElementIndex1Comparator implements Serializable, Comparator<LDPair> {

        public int compare(final LDPair arg0, final LDPair arg1) {
            return arg0.i1 - arg1.i1;
        }
    }
