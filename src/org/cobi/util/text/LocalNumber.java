/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.util.text;

/**
 *
 * @author mxli
 */
public class LocalNumber {

    public static String decimal2XSystem(int d, int xint, int len) {
        int s = 1;
        int n = d;

        while (n >= xint) {
            s++;
            n = n / xint;
        }

        if (len < s) {
            len = s;
        }
        StringBuffer x = new StringBuffer();
        x.setLength(len);
        for (int i = 0; i < len; i++) {
            x.setCharAt(i, '0');
        }
        if (d < xint) {
            x.setCharAt(len - 1, change2X(d));
        } else {
            int c;
            int i = len;
            do {
                i--;
                c = d / xint;
                x.setCharAt(i, change2X(d % xint));//�ж��Ƿ����10���������10����ת��ΪA~F�ĸ�ʽ
                d = c;
            } while (c >= xint);
            x.setCharAt(i - 1, change2X(d));
        }
        return x.toString();
    }

    public static int xSystem2decimal(String x, int xint) {
        int len = x.length();
        int toal = 0;
        for (int i = len - 1; i >= 0; i--) {
            toal += change2Deci(x.charAt(i)) * Math.pow(xint, len - 1 - i);
        }
        return toal;
    }

    public static int xSystem2decimal(StringBuffer x, int xint) {
        int len = x.length();
        int toal = 0;
        for (int i = len - 1; i >= 0; i--) {
            toal += change2Deci(x.charAt(i)) * Math.pow(xint, len - 1 - i);
        }
        return toal;
    }

    //�ж��Ƿ�Ϊ10~15֮�����������������ת��
    private static char change2X(int d) {
        char x = '0';
        switch (d) {
            case 10:
                x = 'A';
                break;
            case 11:
                x = 'B';
                break;
            case 12:
                x = 'C';
                break;
            case 13:
                x = 'D';
                break;
            case 14:
                x = 'E';
                break;
            case 15:
                x = 'F';
                break;
            default:
                x = (char) (d + '0');
                break;
        }
        return x;
    }

    //�ж��Ƿ�Ϊ10~15֮�����������������ת��
    private static int change2Deci(char x) {
        int d = -1;
        switch (x) {
            case 'A':
                d = 10;
                break;
            case 'B':
                d = 11;
                break;
            case 'C':
                d = 12;
                break;
            case 'D':
                d = 13;
                break;
            case 'E':
                d = 14;
                break;
            case 'F':
                d = 15;
                break;
            default:
                d = x - '0';
                break;
        }
        return d;
    }

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        try {
            System.out.println(LocalNumber.decimal2XSystem(8, 2, 9));
            System.out.println(LocalNumber.xSystem2decimal(LocalNumber.decimal2XSystem(154560, 9, 100), 9));
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }
}
