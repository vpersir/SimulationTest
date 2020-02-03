/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.genetsimulator.entity;

import cern.colt.list.IntArrayList;
import cern.colt.matrix.DoubleMatrix2D;
import java.util.Set;

/**
 *
 * @author MX Li
 */
public abstract class LDSparseMatrix {

  //  public abstract float getLDAt(int index1, int index2) throws Exception;

   public abstract double getLDAt(int index1, int index2) throws Exception;

    public abstract DoubleMatrix2D subDenseLDMatrix(IntArrayList indexes) throws Exception;

    public abstract Set<Integer> getAllUniqueIndexes();

    public abstract void releaseLDData();

    public abstract boolean isEmpty();

    public abstract boolean calculateGenotypeCorrelationBlock(IntArrayList indexes) throws Exception;
}
