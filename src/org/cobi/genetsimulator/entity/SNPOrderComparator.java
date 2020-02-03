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
public class SNPOrderComparator implements Comparator {

    public int compare(Object arg0, Object arg1) {
        AnnotSNP obj1 = (AnnotSNP) arg0;
        AnnotSNP obj2 = (AnnotSNP) arg1;
        return (int) (obj1.order - obj2.order);
    }
}
