/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.util.math;

import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author MX Li
 */
public class MultipleSelectionGenerator {

    private int cellNum;
    int[] buffer = null;
    List<int[]> allAnswers = null;
    int needScore;

    public MultipleSelectionGenerator(int cellNum, int needScore) {
        this.cellNum = cellNum;
        this.needScore = needScore;
        buffer = new int[cellNum];
        allAnswers = new ArrayList<int[]>();
    }

    public void f(int[] a, int v, int cellNum) {
        int score = 3;
        if (v != 0) {
            for (int i = 0; i < score; i++) {
                int[] aa = new int[cellNum];
                for (int j = 0; j < cellNum - v; j++) {
                    aa[j] = a[j];
                }
                aa[cellNum - v] = i;
                f(aa, v - 1, cellNum);
            }
        } else {
            for (int j = 0; j < cellNum; j++) {
                System.out.print(a[j]);
                System.out.print(' ');
            }
            System.out.println();
        }
    }

    public void produceAll(int v, int maxScore) {
        if (v < cellNum) {
            for (int i = 0; i < maxScore; i++) {
                buffer[v] = i;
                produceAll(v + 1, maxScore);
            }
        } else {
            for (int j = 0; j < cellNum; j++) {
                System.out.print(buffer[j]);
                System.out.print(' ');
            }
            System.out.println();
        }
    }

    public void producePart(int v, int maxScore) {
        if (v < cellNum) {
            int sum = 0;
            for (int i = 0; i < maxScore; i++) {
                buffer[v] = i;
                sum += i;

                if (sum <= needScore) {
                    producePart(v + 1, maxScore);
                }
            }
        } else {
            int sum = 0;
            for (int j = 0; j < cellNum; j++) {
                sum += buffer[j];
            }
            if (sum == needScore) {
                int[] aa = new int[cellNum];
                System.arraycopy(buffer, 0, aa, 0, cellNum);
                allAnswers.add(aa);
                /*
                for (int j = 0; j < cellNum; j++) {
                System.out.print(buffer[j]);
                System.out.print(' ');
                }
                System.out.println();
                 *
                 */

            }

        }
    }

    public List<int[]> getAllAnswers() {
        return allAnswers;
    }

    public void setAllAnswers(List<int[]> allAnswers) {
        this.allAnswers = allAnswers;
    }

    public static void main(String[] args) {
        int cellNum = 3;
        int maxScore = 5;
        MultipleSelectionGenerator x = new MultipleSelectionGenerator(cellNum, 3);
        //x.produceAll(0, maxScore);
        x.producePart(0, 5);
    }
}
