package org.cobi.eqtlsimulator;

import java.util.Comparator;

public class DoubleArrayListComparatorR implements Comparator<double[]> {

    int index;

    public DoubleArrayListComparatorR(int index) {
        this.index = index;
    }

    @Override
    public int compare(final double[] arg0, final double[] arg1) {
        return Double.compare(arg1[index], arg0[index]);
    }
}
