/*
 * LocalString.java
 *
 * Created on 2008��1��28��, ����4:28
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */
package org.cobi.util.text;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Arrays;
import java.util.Locale;
import org.apache.commons.math.linear.RealMatrix;

/**
 *
 * @author Miaoxin Li
 */
public class LocalString {

    /** Creates a new instance of LocalString */
    public LocalString() {
    }

    /**
     * support Numeric format:<br>
     * "33" "+33" "033.30" "-.33" ".33" " 33." " 000.000 "
     * @param str String
     * @return boolean
     */
    public static boolean isNumeric(String str) {
        int begin = 0;
        boolean once = true;
        if (str == null || str.trim().equals("")) {
            return false;
        }
        str = str.trim();
        if (str.startsWith("+") || str.startsWith("-")) {
            if (str.length() == 1) {
                // "+" "-"
                return false;
            }
            begin = 1;
        }
        //scientific formate like "2.7266453405784747E-4"
        if (str.indexOf("E-") >= 0) {
            str = str.replaceAll("E-", "");
        } else if (str.indexOf("e-") >= 0) {
            str = str.replaceAll("e-", "");
        } else if (str.indexOf("E+") >= 0) {
            str = str.replaceAll("E+", "");
        } else if (str.indexOf("e+") >= 0) {
            str = str.replaceAll("e+", "");
        } else if (str.indexOf("E") >= 0) {
            str = str.replaceAll("E", "");
        } else if (str.indexOf("e") >= 0) {
            str = str.replaceAll("e", "");
        }

        for (int i = begin; i < str.length(); i++) {
            if (!Character.isDigit(str.charAt(i))) {
                if (str.charAt(i) == '.' && once) {
                    // '.' can only once
                    once = false;
                } else {
                    return false;
                }
            }
        }
        if (str.length() == (begin + 1) && !once) {
            // "." "+." "-."
            return false;
        }
        return true;
    }
     public static void print2DRealMatrix(RealMatrix items) {
        int rowNum = items.getRowDimension();
        int colNum = items.getColumnDimension() - 1;
        DecimalFormat df = new DecimalFormat("0.0000000", new DecimalFormatSymbols(Locale.US));
        for (int c = 0; c < rowNum; c++) {
            for (int d = 0; d < colNum; d++) {
                System.out.print(df.format(items.getEntry(c, d)) + "\t");
            }
            System.out.print(df.format(items.getEntry(c, colNum)));
            System.out.println();
        }
    }
    /**
      * support Integer format:<br>
      * "33" "003300" "+33" " -0000 "
      * @param str String
      * @return boolean
      */
    public static boolean isInteger(String str) {
        int begin = 0;
        if (str == null || str.trim().equals("")) {
            return false;
        }
        str = str.trim();
        if (str.startsWith("+") || str.startsWith("-")) {
            if (str.length() == 1) {
                // "+" "-"
                return false;
            }
            begin = 1;
        }
        for (int i = begin; i < str.length(); i++) {
            if (!Character.isDigit(str.charAt(i))) {
                return false;
            }
        }
        return true;
    }


    /**
     * A special function to divide a long String into several sub-strings.
     * Say, "I, love, you,,," will be splited into 6 words here but 4 words
     * by str.split()
     * @param str
     * @param delim
     * @return
     */
    static public String[] splitString(String str, char delim) {
        if (str != null) {
            int len = str.length();
            int blankcell = 0;
            while (str.charAt(--len) == delim) {
                blankcell++;
            }
            String[] cell = str.split(String.valueOf(delim));
            if (blankcell == 0) {
                return cell;
            } else {
                len = cell.length;
                String[] newcell = new String[len + blankcell];
                for (int i = 0; i < len; i++) {
                    newcell[i] = cell[i];
                }
                Arrays.fill(newcell, len, len + blankcell, "");
                return newcell;
            }

        } else {
            return null;
        }
    }

    /**
     * split the String according to regex and ignore the repeat regexs<br>
     * eg. "oooaoodsoooosdaodosoo" regex="o" -> {"a","ds","sda","d","s"}
     * @param src String the source string
     * @param regex String what string splited by
     * @return String[]
     */
    public static String[] splitIgnoreRepeat(String src, String regex) {
        String[] s = src.split(regex);
        String[] d = null;
        int num = 0;
        for (int i = 0; i < s.length; i++) {
            if (!s[i].equals("")) {
                num++;
            }
        }
        d = new String[num];
        num = 0;
        for (int i = 0; i < s.length; i++) {
            if (!s[i].equals("")) {
                d[num++] = s[i];
            }
        }
        return d;
    }
}
