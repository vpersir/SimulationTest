/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.genetsimulator.entity;

import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author mxli
 */
public class GeneExpressSet {

    public String geneSymb;
    public List<double[]> expressSet;

    public GeneExpressSet(String geneSymb) {
        this.geneSymb = geneSymb;
        expressSet = new ArrayList<double[]>();
    }

    public void addExpression(double[] s) {
        expressSet.add(s);
    }
}
