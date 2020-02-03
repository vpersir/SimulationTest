//Source code of HaploView4.0, Copied in April 2008
//Modified by MXLi in May 2008 
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.util.text;

import java.io.File;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.FieldPosition;
import java.text.NumberFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Locale;

/**
 *
 * @author Miaoxin Li
 */
public class Util {

    public static String formatPValue(double pval) {
        DecimalFormat df;
        //java truly sucks for simply restricting the number of sigfigs but still
        //using scientific notation when appropriate
        /*
        if (pval < 0.0001) {
        df = new DecimalFormat("0.000E0", new DecimalFormatSymbols(Locale.US));
        } else {
        df = new DecimalFormat("0.0000000", new DecimalFormatSymbols(Locale.US));
        }
         */
        df = new DecimalFormat("0.00E0", new DecimalFormatSymbols(Locale.US));
        String formattedNumber = df.format(pval, new StringBuffer(), new FieldPosition(NumberFormat.INTEGER_FIELD)).toString();
        return formattedNumber;
    }

    /**
     * Rounds a double and converts it into String.
     *
     * @param value the double value
     * @param afterDecimalPoint the (maximum) number of digits permitted
     * after the decimal point
     * @return the double as a formatted string
     */
    public static String doubleToString(double value, int afterDecimalPoint) {

        StringBuffer stringBuffer;
        double temp;
        int i, dotPosition;
        long precisionValue;

        temp = value * Math.pow(10.0, afterDecimalPoint);
        if (Math.abs(temp) < Long.MAX_VALUE) {
            precisionValue = (temp > 0) ? (long) (temp + 0.5)
                    : -(long) (Math.abs(temp) + 0.5);
            if (precisionValue == 0) {
                stringBuffer = new StringBuffer(String.valueOf(0));
            } else {
                stringBuffer = new StringBuffer(String.valueOf(precisionValue));
            }
            if (afterDecimalPoint == 0) {
                return stringBuffer.toString();
            }
            dotPosition = stringBuffer.length() - afterDecimalPoint;
            while (((precisionValue < 0) && (dotPosition < 1))
                    || (dotPosition < 0)) {
                if (precisionValue < 0) {
                    stringBuffer.insert(1, '0');
                } else {
                    stringBuffer.insert(0, '0');
                }
                dotPosition++;
            }
            stringBuffer.insert(dotPosition, '.');
            if ((precisionValue < 0) && (stringBuffer.charAt(1) == '.')) {
                stringBuffer.insert(1, '0');
            } else if (stringBuffer.charAt(0) == '.') {
                stringBuffer.insert(0, '0');
            }
            int currentPos = stringBuffer.length() - 1;
            while ((currentPos > dotPosition)
                    && (stringBuffer.charAt(currentPos) == '0')) {
                stringBuffer.setCharAt(currentPos--, ' ');
            }
            if (stringBuffer.charAt(currentPos) == '.') {
                stringBuffer.setCharAt(currentPos, ' ');
            }

            return stringBuffer.toString().trim();
        }
        return new String("" + value);
    }

    /**
     * Rounds a double and converts it into a formatted decimal-justified String.
     * Trailing 0's are replaced with spaces.
     *
     * @param value the double value
     * @param width the width of the string
     * @param afterDecimalPoint the number of digits after the decimal point
     * @return the double as a formatted string
     */
    public static String doubleToString(double value, int width,
            int afterDecimalPoint) {

        String tempString = doubleToString(value, afterDecimalPoint);
        char[] result;
        int dotPosition;

        if ((afterDecimalPoint >= width)
                || (tempString.indexOf('E') != -1)) { // Protects sci notation
            return tempString;
        }

        // Initialize result
        result = new char[width];
        for (int i = 0; i < result.length; i++) {
            result[i] = ' ';
        }

        if (afterDecimalPoint > 0) {
            // Get position of decimal point and insert decimal point
            dotPosition = tempString.indexOf('.');
            if (dotPosition == -1) {
                dotPosition = tempString.length();
            } else {
                result[width - afterDecimalPoint - 1] = '.';
            }
        } else {
            dotPosition = tempString.length();
        }


        int offset = width - afterDecimalPoint - dotPosition;
        if (afterDecimalPoint > 0) {
            offset--;
        }

        // Not enough room to decimal align within the supplied width
        if (offset < 0) {
            return tempString;
        }

        // Copy characters before decimal point
        for (int i = 0; i < dotPosition; i++) {
            result[offset + i] = tempString.charAt(i);
        }

        // Copy characters after decimal point
        for (int i = dotPosition + 1; i < tempString.length(); i++) {
            result[offset + i] = tempString.charAt(i);
        }

        return new String(result);
    }

    public static double roundDouble(double d, int places) {
        double factor = Math.pow(10, places);
        return Math.rint(d * factor) / factor;
    }

    public static String getTimeString() {
        /*
        Calendar c = Calendar.getInstance();
        int intYear = c.get(Calendar.YEAR);
        int intMonth = c.get(Calendar.MONTH+1);
        int intDate = c.get(Calendar.DATE);
        int intHour = c.get(Calendar.HOUR_OF_DAY);
        int intMinute = c.get(Calendar.MINUTE);
        int intSecond = c.get(Calendar.SECOND);
        StringBuffer str = new StringBuffer();
        str.append(intYear);
        str.append('_');
        str.append(intMonth);
        str.append('_');
        str.append(intDate);
        str.append('_');
        str.append(intHour);
        str.append('_');
        str.append(intMinute);
        str.append('_');
        str.append(intSecond);
        return str.toString();
         */

        Date d = new Date(System.currentTimeMillis());
        //SimpleDateFormat format = new SimpleDateFormat("yyyy_MM_dd_HH_mm_ss");
        SimpleDateFormat format = new SimpleDateFormat("MM_dd_HH_mm_ss");
        return format.format(d);
    }

    public static boolean makeStorageLoc(String filePath) throws Exception {
        if (filePath == null) {
            return false;
        }

        File nwFile = new File(filePath);
        if (!nwFile.exists()) {
            if (!nwFile.mkdirs()) {
                return false;
            }
        }
        return true;
    }

   
}

 