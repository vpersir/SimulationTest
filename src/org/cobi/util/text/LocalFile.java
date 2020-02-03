/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.util.text;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.StringTokenizer;
import java.util.HashSet;

/**
 *
 * @author MX Li
 */
public class LocalFile {

    /**
     * retrieve data from a text file whith limited rows
     * @param FileName
     * @param indices
     * @param arry
     * @throws java.lang.Exception
     * @return
     */
    static public boolean retrieveData(String fileName, List<String[]> arry,
            int limitedRowNumber, String delimi, String startLabel, boolean useTokenizer) throws Exception {
        BufferedReader br = new BufferedReader(new FileReader(fileName));
        String line = "";
        String[] row = null;
        int lineNumber = 0;
        String delmilit = "\t\" \"\n,";
        if (delimi != null) {
            delmilit = delimi;
            //usually some files donot start with data but with breif annoation, so we need filter the latter.
        }
        if (startLabel != null) {
            while ((line = br.readLine()) != null) {
                if (line.trim().startsWith(startLabel)) {
                    break;
                }
            }
        }


        if (useTokenizer) {
            int colNum = -1;
            int i;
            StringBuilder tmpStr = new StringBuilder();
            do {
                if (line.trim().length() == 0) {
                    continue;
                }
                StringTokenizer tokenizer = new StringTokenizer(line, delmilit);
                if (colNum < 0) {
                    colNum = tokenizer.countTokens();
                }
                row = new String[colNum];
                for (i = 0; i < colNum; i++) {
                    //sometimes tokenizer.nextToken() can not release memory
                    row[i] = tmpStr.append(tokenizer.nextToken().trim()).toString();
                    tmpStr.delete(0, tmpStr.length());
                }
                arry.add(row);

                lineNumber++;
                if (lineNumber > limitedRowNumber) {
                    break;
                }
            } while ((line = br.readLine()) != null);
        } else {
            if (delmilit.equals("\t\" \"\n")) {
                delmilit = "[" + delmilit + "]";
            }
            do {
                if (line.trim().length() == 0) {
                    continue;
                }
                arry.add(line.split(delmilit, -1));
                lineNumber++;
                if (lineNumber > limitedRowNumber) {
                    break;
                }
            } while ((line = br.readLine()) != null);
        }
        br.close();
        return true;
    }

    /**
     * simply retrieve data from a  file
     * @param FileName
     * @param arry
     * @throws java.lang.Exception
     * @return
     */
    static public boolean retrieveData(String fileName, List<String[]> arry, String delimiter) throws Exception {
        BufferedReader br = new BufferedReader(new FileReader(fileName));
        String line = null;
        String delmilit = "\t\" \"\n";
        if (delimiter != null) {
            delmilit = delimiter;        //usually some files donot start with data but with breif annoation, so we need filter the latter.

        }
        int colNum = -1;
        String[] row = null;
        StringBuilder tmpStr = new StringBuilder();
        while ((line = br.readLine()) != null) {
            if (line.trim().length() == 0) {
                continue;
            }
            StringTokenizer tokenizer = new StringTokenizer(line, delmilit);
            if (colNum < 0) {
                colNum = tokenizer.countTokens();
            }
            row = new String[colNum];
            for (int i = 0; i < colNum; i++) {
                //sometimes tokenizer.nextToken() can not release memory
                row[i] = tmpStr.append(tokenizer.nextToken().trim()).toString();
                tmpStr.delete(0, tmpStr.length());
            }
            arry.add(row);
        }
        br.close();
        return true;
    }

    /**
     * simply retrieve data from a  file
     * @param FileName
     * @param arry
     * @throws java.lang.Exception
     * @return
     */
    static public boolean retrieveData(String fileName, StringBuilder tmpBf) throws Exception {
        BufferedReader br = new BufferedReader(new FileReader(fileName));
        String line = null;
        while ((line = br.readLine()) != null) {
            line = line.trim();
            tmpBf.append(line);
            tmpBf.append('\n');
        }
        br.close();
        return true;
    }

    /**
     * simply retrieve data from a  file
     * @param FileName
     * @param arry
     * @throws java.lang.Exception
     * @return
     */
    static public boolean retrieveData(String fileName, HashSet<String> arry) throws Exception {
        BufferedReader br = new BufferedReader(new FileReader(fileName));
        String line = null;

        while ((line = br.readLine()) != null) {
            line = line.trim();
            if (line.length() > 1) {
                arry.add(line);
            }
        }
        br.close();
        return true;
    }

    /**
     * simply retrieve data from a  file
     * @param FileName
     * @param arry
     * @throws java.lang.Exception
     * @return
     */
    static public boolean retrieveData(String fileName, HashSet<String> arry, int index) throws Exception {
        BufferedReader br = new BufferedReader(new FileReader(fileName));
        String line = null;
        int i = 0;
        while ((line = br.readLine()) != null) {
            StringTokenizer tokenizer = new StringTokenizer(line);
            for (i = 0; i < index; i++) {
                tokenizer.nextToken();
            }
            arry.add(tokenizer.nextToken()); 
        }
        br.close();
        return true;
    }

    /**
     * retrieve data from a text file it is based on split,
     * @param FileName
     * @param indices
     * @param arry
     * @throws java.lang.Exception
     * @return
     */
    static public boolean retrieveData(String fileName, List<String[]> arry, int[] orgIndices,
            HashSet<String> refList, int refIndex, String delimiter) throws Exception {
        BufferedReader br = new BufferedReader(new FileReader(fileName));
        String line = null;
        String[] cells = null;
        String[] row = null;
        int selectedColNum = orgIndices.length;

        int i, pos;

        while ((line = br.readLine()) != null) {
            line = line.trim();
            if (line.trim().length() == 0) {
                continue;
            }
            cells = line.split(delimiter, -1);
            if (refList.contains(cells[refIndex])) {
                row = new String[selectedColNum];
                for (i = 0; i < selectedColNum; i++) {
                    row[i] = cells[orgIndices[i]];
                }
                arry.add(row);
            }
        }
        br.close();
        return true;
    }

    /**
     * retrieve data from a text file it is based on tokenizor, but the order in the indices will
     * be changed to the order of the files. it is a consideration of speed.
     * @param FileName
     * @param indices
     * @param arry
     * @throws java.lang.Exception
     * @return
     */
    static public boolean retrieveData(String fileName, List<String[]> arry, int[] orgIndices,
            String delimiter) throws Exception {
        BufferedReader br = new BufferedReader(new FileReader(fileName));
        String line = null;
        String[] cells = null;
        String[] row = null;
        int selectedColNum = orgIndices.length;
        int i;

        while ((line = br.readLine()) != null) {
            //line = line.trim();
            if (line.trim().length() == 0) {
                continue;
            }
            cells = line.split(delimiter, -1);
            row = new String[selectedColNum];
            for (i = 0; i < selectedColNum; i++) {
                row[i] = cells[orgIndices[i]];
            }
            arry.add(row);
        }
        br.close();
        return true;
    }

    /**
     * retrieve data from a text file it is based on split,
     * @param FileName
     * @param indices
     * @param arry
     * @throws java.lang.Exception
     * @return
     */
    static public boolean retrieveData(String fileName, ArrayList<String[]> arry, int[] orgIndices,
            String[] refList, int refIndex, String delimiter) throws Exception {
        BufferedReader br = new BufferedReader(new FileReader(fileName));
        String line = null;
        String[] cells = null;
        String[] row = null;
        int selectedColNum = orgIndices.length;
        int i, pos;
        Arrays.sort(refList);
        while ((line = br.readLine()) != null) {
            line = line.trim();
            if (line.trim().length() == 0) {
                continue;
            }
            cells = line.split(delimiter, -1);
            pos = Arrays.binarySearch(refList, cells[refIndex]);
            if (pos >= 0) {
                row = new String[selectedColNum];
                for (i = 0; i < selectedColNum; i++) {
                    row[i] = cells[orgIndices[i]];
                }
                arry.add(row);
            }
        }
        br.close();
        return true;
    }

    /**
     * write data to a text file
     * @param FileName
     * @param indices
     * @param arry
     * @throws java.lang.Exception
     * @return
     */
    static public boolean writeObject2Text(String fileName, List<Object[]> arry, String delmilit) throws Exception {
        BufferedWriter bw = new BufferedWriter(new FileWriter(fileName));
        Object[] linecells = null;
        int linenumber = arry.size();
        int cols = 0;
        for (int i = 0; i < linenumber; i++) {
            linecells = arry.get(i);
            cols = linecells.length - 1;
            for (int j = 0; j < cols; j++) {
                if (linecells[j] == null) {
                    bw.write(" ");
                } else {
                    bw.write(linecells[j].toString());
                }
                bw.write("\t");
            }
            if (linecells[cols] == null) {
                bw.write(" ");
            } else {
                bw.write(linecells[cols].toString());
            }
            bw.write("\n");

        }
        bw.flush();
        bw.close();
        return true;
    }

    /**
     * write data to a text file
     * @param FileName
     * @param indices
     * @param arry
     * @throws java.lang.Exception
     * @return
     */
    static public boolean writeData(String fileName, List<String[]> arry, String delmilit, boolean append) throws Exception {
        BufferedWriter bw = new BufferedWriter(new FileWriter(fileName, append));
        String[] linecells = null;
        int linenumber = arry.size();
        if (linenumber == 0) {
            return false;
        }
        int cols = arry.get(0).length - 1;
        for (int i = 0; i < linenumber; i++) {
            linecells = arry.get(i);

            for (int j = 0; j < cols; j++) {
                if (linecells[j] == null) {
                    bw.write("-");
                } else {
                    bw.write(linecells[j]);
                }
                bw.write(delmilit);
            }

            if (linecells[cols] == null) {
                bw.write("-");
            } else {
                bw.write(linecells[cols]);
            }

            bw.write("\n");
        }
        bw.close();
        return true;
    }
}
