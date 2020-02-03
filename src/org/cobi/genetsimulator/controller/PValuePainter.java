/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.genetsimulator.controller;

import cern.colt.list.DoubleArrayList;
import cern.jet.stat.Descriptive;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GradientPaint;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.RenderingHints;
import java.awt.Stroke;
import java.awt.geom.Line2D;
import java.awt.geom.Path2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.Iterator;
import java.util.List;
import javax.imageio.IIOImage;
import javax.imageio.ImageIO;
import javax.imageio.ImageTypeSpecifier;
import javax.imageio.ImageWriteParam;
import javax.imageio.ImageWriter;
import javax.imageio.stream.ImageOutputStream;
import org.cobi.util.text.Util;

/**
 *
 * @author mxli
 */
public class PValuePainter { //unit pixel

    private int offsetLeft = 70;
    private int offsetTop = 10;
    private int offsetBottom = 40;
    private int offsetRight = 20;
    private Rectangle plotScope;
    private Dimension canvasScope = new Dimension(650, 400);
    private boolean hasTitleFlag = false;
    private Color plotBackgroundColor = Color.WHITE;
    private Color axesColor = Color.DARK_GRAY;
    private Color nonplotBackgroundColor = Color.LIGHT_GRAY;
    private Font titleFont = new Font("Dialog", Font.PLAIN, 16);

    /**
     *
     * @param width
     * @param height
     */
    public PValuePainter(int width, int height) {
        canvasScope.width = width;
        canvasScope.height = height;
        calculatePlotArea();
    }

    /**
     *
     * @param valueList
     * @param title
     * @param outputPath
     * @param pValueTolerationLevle
     * @return
     * @throws Exception
     */
    public int drawQQPlot(final DoubleArrayList valueList, String title, String outputPath, double pValueTolerationLevle) throws Exception {
        if (valueList == null || valueList.size() == 0) {
            System.err.println("Null p-value list");
            return -1;
        }
        DoubleArrayList tmpValueList = valueList.copy();
        BufferedImage image = new BufferedImage(this.canvasScope.width, this.canvasScope.height, BufferedImage.TYPE_INT_RGB);
        Graphics2D g2d = image.createGraphics();
        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        paintPlotArea(g2d);
        drawAxes(g2d, title, "Expected (-log(10)P)", "Observed (-log(10)P)");

        int dataSize = tmpValueList.size();
        tmpValueList.quickSort();
        //max and min value in vertical
        double vMax = -Math.log10(tmpValueList.getQuick(0));
        double vMin = -Math.log10(tmpValueList.getQuick(dataSize - 1));
        double hMin = 0;
        double hMax = -Math.log10((1.0 / (dataSize + 1)));

        //sometimes the p values are too large
        if (vMax > -Math.log10(pValueTolerationLevle)) {
            vMax = -Math.log10(pValueTolerationLevle);
        }

        CoordinateTransformer transformer = new CoordinateTransformer();
        transformer.setupBasicScope(plotScope.x, plotScope.x + plotScope.width,
                plotScope.y, plotScope.y + plotScope.height, hMin, hMax, vMin, vMax);
        drawAxesScale(g2d, transformer);

        //draw the standard line
        Point2D point1 = transformer.data2ScreenPoint(0, 0);
        Point2D point2 = transformer.data2ScreenPoint(Math.max(hMax, vMax), Math.max(hMax, vMax));
        Line2D.Double zz = new Line2D.Double(point1, point2);
        Stroke oldStroke = g2d.getStroke();
        Color oldColor = g2d.getColor();
        g2d.setColor(Color.GREEN);
        g2d.setStroke(new BasicStroke(2f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
        g2d.draw(zz);
        g2d.setStroke(oldStroke);
        g2d.setColor(oldColor);

        double rectangleSize = 1;
        double maxDiff = 0.05;
        int maxDiffIndex = -1;
        double tmpDouble2;
        double tmpDouble1;
        double triangleSide = 4;
        double halfTriangleSide = triangleSide / 2;
        double triangleHeight = Math.cos(Math.PI / 3) * triangleSide;
        int dataSizeMore = dataSize + 1;
        for (int i = 0; i < dataSize; i++) {
            tmpDouble1 = -Math.log10((double) (i + 1) / (dataSizeMore));
            tmpDouble2 = -Math.log10(tmpValueList.getQuick(i));
            if (tmpDouble2 <= vMax) {
                point1 = transformer.data2ScreenPoint(tmpDouble1, tmpDouble2);
                Rectangle2D.Double f = new Rectangle2D.Double(point1.getX() - rectangleSize / 2,
                        point1.getY() - rectangleSize / 2, rectangleSize, rectangleSize);
                g2d.draw(f);
            } else {
                tmpDouble2 = vMax;
                point1 = transformer.data2ScreenPoint(tmpDouble1, tmpDouble2);
                double x = point1.getX();
                double y = point1.getY();
                double[] xs = {x, x + halfTriangleSide, x - halfTriangleSide};
                double[] ys = {y, y - triangleHeight, y - triangleHeight};
                Path2D.Double triangle = new Path2D.Double();
                triangle.moveTo(xs[0], ys[0]);
                for (int j = 1; j < xs.length; j++) {
                    triangle.lineTo(xs[j], ys[j]);
                }
                triangle.lineTo(xs[0], ys[0]);

                g2d.setColor(Color.red);
                g2d.draw(triangle);
                g2d.setColor(oldColor);
            }
            /*
            if (maxDiffIndex < 0 && (Math.abs(expectedList.getQuick(i) - tmpValueList.getQuick(i)) > maxDiff)) {
            maxDiffIndex = i;
            }
             *
             */
        }

        /*
        float[] dashes = {3.f};
        g2d.setColor(Color.RED);
        g2d.setStroke(new BasicStroke(0.5f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND, 10, dashes, 0));

        //draw large difference
        point1 = transformer.data2ScreenPoint(expectedList.getQuick(maxDiffIndex), tmpValueList.getQuick(maxDiffIndex));
        point2 = transformer.data2ScreenPoint(hMax, vMax);
        double vRectangleSize = point1.getY() - point2.getY();
        double hRectangleSize = point2.getX() - point1.getX();
        //zz = new Line2D.Double(point1, point2);
        // g2d.draw(zz);
        Rectangle2D.Double f = new Rectangle2D.Double(point1.getX(),
        point1.getY() - vRectangleSize, hRectangleSize, vRectangleSize);
        g2d.draw(f);

        g2d.setStroke(oldStroke);
        g2d.setColor(oldColor);
         */

        g2d.dispose();
        //outputJPEGFile(image, outputPath);
        outputPNGFile(image, outputPath);
        return maxDiffIndex;
    }

    /**
     *
     * @param valueLists
     * @param legends
     * @param title
     * @param outputPath
     * @param pValueTolerationLevle
     * @return
     * @throws Exception
     */
    public int drawMultipleQQPlot(final List<DoubleArrayList> valueLists, List<String> legends,
            String title, String outputPath, double pValueTolerationLevle) throws Exception {
        if (valueLists == null || valueLists.isEmpty()) {
            System.err.println("Null p-value list");
            return -1;
        }

        int listNum = valueLists.size();
        for (int i = listNum - 1; i >= 0; i--) {
            if (valueLists.get(i) == null || valueLists.get(i).size() == 0) {
                System.err.println("Null p-value list");
                valueLists.remove(i);
            }
        }

        if (valueLists == null || valueLists.isEmpty()) {
            System.err.println("Null p-value list");
            return -1;
        }

        int dataSize = valueLists.get(0).size();

        //max and min value in vertical
        double vMax = (valueLists.get(0).getQuick(0));
        double vMin = vMax;
        double hMin = 0;
        double hMax = (1.0 / dataSize);
        double tmpDouble;

        for (int i = 0; i < listNum; i++) {
            DoubleArrayList tmpValueList = valueLists.get(i);
            tmpValueList.quickSort();
            dataSize = tmpValueList.size();
            for (int j = 0; j < dataSize; j++) {
                if (vMax > tmpValueList.getQuick(j)) {
                    vMax = tmpValueList.getQuick(j);
                } else if (vMin < tmpValueList.getQuick(j)) {
                    vMin = tmpValueList.getQuick(j);
                }
            }
            if (hMax > (1.0 / (dataSize + 1))) {
                hMax = (1.0 / (dataSize + 1));
            }
        }

        //convert them into the data to be used
        vMin = -Math.log10(vMin);
        vMax = -Math.log10(vMax);
        hMax = -Math.log10(hMax);

        //sometimes the p values are too large
        if (vMax > -Math.log10(pValueTolerationLevle)) {
            vMax = -Math.log10(pValueTolerationLevle);
        }

        BufferedImage image = new BufferedImage(canvasScope.width, canvasScope.height, BufferedImage.TYPE_INT_RGB);
        Graphics2D g2d = image.createGraphics();

        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        paintPlotArea(g2d);
        drawAxes(g2d, title, "Expected (-log(10)P)", "Observed (-log(10)P)");
        CoordinateTransformer transformer = new CoordinateTransformer();
        transformer.setupBasicScope(plotScope.x, plotScope.x + plotScope.width,
                plotScope.y, plotScope.y + plotScope.height, hMin, hMax, vMin, vMax);
        drawAxesScale(g2d, transformer);
        //draw the standard line
        Point2D point1 = transformer.data2ScreenPoint(0, 0);
        Point2D point2 = transformer.data2ScreenPoint(Math.max(hMax, vMax), Math.max(hMax, vMax));
        Line2D.Double zz = new Line2D.Double(point1, point2);
        Stroke oldStroke = g2d.getStroke();
        Color oldColor = g2d.getColor();
        g2d.setColor(Color.GREEN);
        g2d.setStroke(new BasicStroke(2f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
        g2d.draw(zz);

        double rectangleSize = 1;
        double triangleSide = 4;
        double halfTriangleSide = triangleSide / 2;
        double triangleHeight = Math.cos(Math.PI / 3) * triangleSide;

        double maxDiff = 0.05;
        int maxDiffIndex = -1;
        Color[] colors = {Color.BLUE, Color.RED, Color.ORANGE, Color.MAGENTA, Color.BLACK, Color.YELLOW};

        //Color[] colors = {new Color(00, 00, 0xff), new Color(0xff, 00, 00), new Color(00, 0x88, 0xff),    new Color(0xff, 0x88, 00), new Color(0x88, 00, 0xff), new Color(00, 00, 00)};
        int str_hei = g2d.getFontMetrics().getAscent();
        double vRectangleSize = 10;
        double hRectangleSize = 10;
        double tmpDouble1;
        double tmpDouble2;
        g2d.setStroke(oldStroke);

        for (int i = listNum - 1; i >= 0; i--) {
            DoubleArrayList tmpValueList = valueLists.get(i);
            dataSize = tmpValueList.size();
            g2d.setColor(colors[i]);
            int dataSizeMore = dataSize + 1;
            for (int j = 0; j < dataSize; j++) {
                tmpDouble1 = -Math.log10((double) (j + 1) / (dataSizeMore));
                tmpDouble2 = -Math.log10(tmpValueList.getQuick(j));
                if (maxDiffIndex < 0 && (Math.abs(tmpDouble1 - tmpValueList.getQuick(j)) > maxDiff)) {
                    maxDiffIndex = j;
                }

                if (tmpDouble2 <= vMax) {
                    point1 = transformer.data2ScreenPoint(tmpDouble1, tmpDouble2);
                    Rectangle2D.Double f = new Rectangle2D.Double(point1.getX() - rectangleSize / 2,
                            point1.getY() - rectangleSize / 2, rectangleSize, rectangleSize);
                    g2d.draw(f);
                } else {
                    tmpDouble2 = vMax;
                    point1 = transformer.data2ScreenPoint(tmpDouble1, tmpDouble2);
                    double x = point1.getX();
                    double y = point1.getY();
                    double[] xs = {x, x + halfTriangleSide, x - halfTriangleSide};
                    double[] ys = {y, y - triangleHeight, y - triangleHeight};
                    Path2D.Double triangle = new Path2D.Double();
                    triangle.moveTo(xs[0], ys[0]);
                    for (int k = 1; k < xs.length; k++) {
                        triangle.lineTo(xs[k], ys[k]);
                    }
                    triangle.lineTo(xs[0], ys[0]);
                    g2d.draw(triangle);
                }

            }
            float[] dashes = {3.f};
            int str_len = g2d.getFontMetrics().stringWidth(legends.get(i));


            //g2d.setStroke(new BasicStroke(0.5f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND, 10, dashes, 0));
            g2d.drawString(legends.get(i), (float) (plotScope.x + plotScope.getWidth() / 4) - str_len / 2,
                    (float) (plotScope.getY() + plotScope.getHeight() / 4) + str_hei * i);
            Rectangle2D.Double f = new Rectangle2D.Double(plotScope.x + plotScope.getWidth() / 4 + str_len / 2 + hRectangleSize,
                    plotScope.getY() + plotScope.getHeight() / 4 - vRectangleSize + str_hei * i, hRectangleSize, vRectangleSize);

            g2d.fill(f);
            g2d.draw(f);

            /*
            g2d.setColor(Color.RED);
            g2d.setStroke(new BasicStroke(0.5f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND, 10, dashes, 0));
            //draw large difference
            point1 = transformer.data2ScreenPoint(-Math.log10((double) (maxDiffIndex + 1) / (dataSize + 1)), tmpValueList.getQuick(maxDiffIndex));
            point2 = transformer.data2ScreenPoint(hMax, vMax);
            double vRectangleSize = point1.getY() - point2.getY();
            double hRectangleSize = point2.getX() - point1.getX();
            //zz = new Line2D.Double(point1, point2);
            // g2d.draw(zz);
            Rectangle2D.Double f = new Rectangle2D.Double(point1.getX(),
            point1.getY() - vRectangleSize, hRectangleSize, vRectangleSize);
            g2d.draw(f);
             */
            g2d.setStroke(oldStroke);
            tmpValueList = null;
        }

        g2d.setColor(oldColor);

        g2d.dispose();
        //outputJPEGFile(image, outputPath);
        outputPNGFile(image, outputPath);
        return maxDiffIndex;
    }

    
    private void calculatePlotArea() {
        try {
            if (isHasTitleFlag()) {
                plotScope = new Rectangle(offsetLeft, offsetTop + 60,
                        canvasScope.width - offsetRight - offsetLeft, canvasScope.height - offsetBottom - offsetTop - 60);
            } else {
                plotScope = new Rectangle(offsetLeft, offsetTop,
                        canvasScope.width - offsetRight - offsetLeft, canvasScope.height - offsetBottom - offsetTop);
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     *
     * @return
     */
    public boolean isHasTitleFlag() {
        return hasTitleFlag;
    }

    /**
     *
     * @param hasTitleFlag
     */
    public void setHasTitleFlag(boolean hasTitleFlag) {
        this.hasTitleFlag = hasTitleFlag;
    }

    /**
     *
     * @param image
     * @param outputPath
     * @throws Exception
     */
    public static void outputJPEGFile(BufferedImage image, String outputPath) throws Exception {
        ImageWriter writer = null;
        ImageTypeSpecifier type = ImageTypeSpecifier.createFromRenderedImage(image);
        Iterator iter = ImageIO.getImageWriters(type, "jpg");
        if (iter.hasNext()) {
            writer = (ImageWriter) iter.next();
        }
        if (writer == null) {
            return;
        }
        IIOImage iioImage = new IIOImage(image, null, null);
        ImageWriteParam param = writer.getDefaultWriteParam();

        param.setCompressionMode(ImageWriteParam.MODE_EXPLICIT);
        param.setCompressionQuality((float) 1.0);
        ImageOutputStream outputStream = ImageIO.createImageOutputStream(new File(outputPath));
        writer.setOutput(outputStream);
        writer.write(null, iioImage, param);
        outputStream.flush();
        outputStream.close();
        writer.dispose();
    }

    /**
     *
     * @param image
     * @param outputPath
     * @throws Exception
     */
    public static void outputPNGFile(BufferedImage image, String outputPath) throws Exception {
        ImageWriter writer = null;
        ImageTypeSpecifier type = ImageTypeSpecifier.createFromRenderedImage(image);
        Iterator iter = ImageIO.getImageWriters(type, "png");
        if (iter.hasNext()) {
            writer = (ImageWriter) iter.next();
        }
        if (writer == null) {
            return;
        }
        IIOImage iioImage = new IIOImage(image, null, null);
        ImageWriteParam param = writer.getDefaultWriteParam();
        ImageOutputStream outputStream = ImageIO.createImageOutputStream(new File(outputPath));
        writer.setOutput(outputStream);
        writer.write(null, iioImage, param);
        outputStream.flush();
        outputStream.close();
        writer.dispose();
    }

    private void paintPlotArea(Graphics2D g2d) {
        g2d.clipRect(0, 0, canvasScope.width, canvasScope.height);
        g2d.setColor(plotBackgroundColor);
        Rectangle2D rc = new Rectangle2D.Double((double) 0,
                (double) 0, (double) canvasScope.width,
                (double) canvasScope.height);
        g2d.fill(rc);
        GradientPaint pat = new GradientPaint(0f, 0f, Color.WHITE, 100f, 45f, Color.GREEN);
        g2d.setPaint(pat);
        /*
        //		fill background
        g2d.setColor(plotBackgroundColor);
        rc = new Rectangle2D.Double((double) plotScope.getX(),
        (double) plotScope.getY(), (double) plotScope.getWidth(),
        (double) plotScope.getHeight());
        g2d.fill(rc);
         *
         */
    }

    /**

     * @param g2
     */
    private void drawAxes(Graphics2D g2, String title, String XLabel, String yLabel) {
        g2.setColor(axesColor.darker());
        //vertical axis
        Line2D.Double zz = new Line2D.Double((double) plotScope.getX(),
                (double) (plotScope.getY()), (double) plotScope.getX(), (double) (plotScope.getY() + plotScope.getHeight()));

        g2.draw(zz);
        //draw arrow
        // g2.drawLine(plotScope.x, plotScope.y, plotScope.x - 3, plotScope.y + 7);
        // g2.drawLine(plotScope.x, plotScope.y, plotScope.x + 3, plotScope.y + 7);

        int str_len = g2.getFontMetrics().stringWidth(yLabel);
        int str_hei = g2.getFontMetrics().getAscent();
        if (str_len > offsetLeft) {
            str_len = offsetLeft;
        }
        g2.drawString(yLabel, plotScope.x - str_len, plotScope.y + str_hei + 5);

        //horizontal axis
        zz = new Line2D.Double((double) plotScope.getX(),
                (double) (plotScope.getY() + plotScope.getHeight()),
                (double) (plotScope.getX() + plotScope.getWidth()),
                (double) (plotScope.getY() + plotScope.getHeight()));
        g2.draw(zz);
        // g2.drawLine(plotScope.x + plotScope.width, plotScope.y + plotScope.height, plotScope.x + plotScope.width - 7, plotScope.y + plotScope.height - 3);
        // g2.drawLine(plotScope.x + plotScope.width, plotScope.y + plotScope.height, plotScope.x + plotScope.width - 7, plotScope.y + plotScope.height + 3);


        str_len = g2.getFontMetrics().stringWidth(XLabel);
        g2.drawString(XLabel, plotScope.x + plotScope.width / 2 - str_len / 2,
                (float) (plotScope.getY() + plotScope.getHeight() + offsetBottom * 2 / 3));


        Font oldFront = g2.getFont();
        g2.setFont(titleFont);
        str_len = g2.getFontMetrics().stringWidth(title);
        g2.drawString(title, plotScope.x + plotScope.width / 2 - str_len / 2,
                (float) (plotScope.getY() + plotScope.getHeight() / 6));
        g2.setFont(oldFront);
    }

    private void drawAxesSciScale(Graphics2D g2, CoordinateTransformer transformer) {
        g2.setColor(axesColor.darker());
        String str;
        int scaleLen = 3;
        int verticalScale = 10;
        double vMin = transformer.getDataVerticalMin();
        double vMax = transformer.getDataVerticalMax();
        double hMin = transformer.getDataHorizontalMin();
        double hMax = transformer.getDataHorizontalMax();
        double gridLen = (vMax - vMin) / verticalScale;
        int str_hei = g2.getFontMetrics().getAscent();
        for (int i = 1; i <= verticalScale; i++) {
            Point2D point1 = transformer.data2ScreenPoint(hMin, vMin + i * gridLen);
            Point2D point2 = new Point2D.Double(point1.getX() + scaleLen, point1.getY());
            Line2D.Double zz = new Line2D.Double(point1, point2);
            g2.draw(zz);
            str = Util.formatPValue(vMin + i * gridLen);

            int str_len = g2.getFontMetrics().stringWidth(str);
            g2.drawString(str, (float) point1.getX() - str_len,
                    (float) (point1.getY() + str_hei / 2));
        }

        int horizontalScale = 10;
        gridLen = (hMax - hMin) / horizontalScale;
        for (int i = 1; i <= horizontalScale; i++) {
            Point2D point1 = transformer.data2ScreenPoint(hMin + i * gridLen, vMin);
            Point2D point2 = new Point2D.Double(point1.getX(), point1.getY() - scaleLen);
            Line2D.Double zz = new Line2D.Double(point1, point2);
            g2.draw(zz);
            str = Util.formatPValue(hMin + i * gridLen);
            int str_len = g2.getFontMetrics().stringWidth(str);
            g2.drawString(str, (float) point1.getX() - str_len / 2,
                    (float) (point1.getY() + str_hei));
        }
    }

    private void drawAxesScale(Graphics2D g2, CoordinateTransformer transformer) {
        g2.setColor(axesColor.darker());
        String str;
        int scaleLen = 3;
        int verticalScale = 10;
        double vMin = transformer.getDataVerticalMin();
        double vMax = transformer.getDataVerticalMax();
        double hMin = transformer.getDataHorizontalMin();
        double hMax = transformer.getDataHorizontalMax();
        double gridLen = (vMax - vMin) / verticalScale;
        int str_hei = g2.getFontMetrics().getAscent();
        for (int i = 1; i <= verticalScale; i++) {
            Point2D point1 = transformer.data2ScreenPoint(hMin, vMin + i * gridLen);
            Point2D point2 = new Point2D.Double(point1.getX() + scaleLen, point1.getY());
            Line2D.Double zz = new Line2D.Double(point1, point2);
            g2.draw(zz);
            str = Util.doubleToString(vMin + i * gridLen, 5, 2);

            int str_len = g2.getFontMetrics().stringWidth(str);
            g2.drawString(str, (float) point1.getX() - str_len,
                    (float) (point1.getY() + str_hei / 2));
        }

        int horizontalScale = 10;
        gridLen = (hMax - hMin) / horizontalScale;
        for (int i = 1; i <= horizontalScale; i++) {
            Point2D point1 = transformer.data2ScreenPoint(hMin + i * gridLen, vMin);
            Point2D point2 = new Point2D.Double(point1.getX(), point1.getY() - scaleLen);
            Line2D.Double zz = new Line2D.Double(point1, point2);
            g2.draw(zz);
            str = Util.doubleToString(hMin + i * gridLen, 5, 2);
            int str_len = g2.getFontMetrics().stringWidth(str);
            g2.drawString(str, (float) point1.getX() - str_len / 2,
                    (float) (point1.getY() + str_hei));
        }
    }

    /**
     *
     */
    public class CoordinateTransformer {

        private int screenHorizontalMin, screenVerticalMin, screenHorizontalMax, screenVerticalMax;
        private double dataHorizontalMin, dataVerticalMin, dataHorizontalMax, dataVerticalMax;
        double horizontalResolution, verticalResolution;

        /**
         *
         * @return
         */
        public double getDataHorizontalMax() {
            return dataHorizontalMax;
        }

        /**
         *
         * @param dataHorizontalMax
         */
        public void setDataHorizontalMax(double dataHorizontalMax) {
            this.dataHorizontalMax = dataHorizontalMax;
        }

        /**
         *
         * @return
         */
        public double getDataHorizontalMin() {
            return dataHorizontalMin;
        }

        /**
         *
         * @param dataHorizontalMin
         */
        public void setDataHorizontalMin(double dataHorizontalMin) {
            this.dataHorizontalMin = dataHorizontalMin;
        }

        /**
         *
         * @return
         */
        public double getDataVerticalMax() {
            return dataVerticalMax;
        }

        /**
         *
         * @param dataVerticalMax
         */
        public void setDataVerticalMax(double dataVerticalMax) {
            this.dataVerticalMax = dataVerticalMax;
        }

        /**
         *
         * @return
         */
        public double getDataVerticalMin() {
            return dataVerticalMin;
        }

        /**
         *
         * @param dataVerticalMin
         */
        public void setDataVerticalMin(double dataVerticalMin) {
            this.dataVerticalMin = dataVerticalMin;
        }

        /**
         *
         * @return
         */
        public double getHorizontalResolution() {
            return horizontalResolution;
        }

        /**
         *
         * @param horizontalResolution
         */
        public void setHorizontalResolution(double horizontalResolution) {
            this.horizontalResolution = horizontalResolution;
        }

        /**
         *
         * @return
         */
        public int getScreenHorizontalMax() {
            return screenHorizontalMax;
        }

        /**
         *
         * @param screenHorizontalMax
         */
        public void setScreenHorizontalMax(int screenHorizontalMax) {
            this.screenHorizontalMax = screenHorizontalMax;
        }

        /**
         *
         * @return
         */
        public int getScreenHorizontalMin() {
            return screenHorizontalMin;
        }

        /**
         *
         * @param screenHorizontalMin
         */
        public void setScreenHorizontalMin(int screenHorizontalMin) {
            this.screenHorizontalMin = screenHorizontalMin;
        }

        /**
         *
         * @return
         */
        public int getScreenVerticalMax() {
            return screenVerticalMax;
        }

        /**
         *
         * @param screenVerticalMax
         */
        public void setScreenVerticalMax(int screenVerticalMax) {
            this.screenVerticalMax = screenVerticalMax;
        }

        /**
         *
         * @return
         */
        public int getScreenVerticalMin() {
            return screenVerticalMin;
        }

        /**
         *
         * @param screenVerticalMin
         */
        public void setScreenVerticalMin(int screenVerticalMin) {
            this.screenVerticalMin = screenVerticalMin;
        }

        /**
         *
         * @return
         */
        public double getVerticalResolution() {
            return verticalResolution;
        }

        /**
         *
         * @param verticalResolution
         */
        public void setVerticalResolution(double verticalResolution) {
            this.verticalResolution = verticalResolution;
        }

        /**
         *
         * @param screenHorizontalMin
         * @param screenHorizontalMax
         * @param screenVerticalMin
         * @param screenVerticalMax
         * @param dataHorizontalMin
         * @param dataHorizontalMax
         * @param dataVerticalMin
         * @param dataVerticalMax
         */
        public void setupBasicScope(int screenHorizontalMin, int screenHorizontalMax, int screenVerticalMin, int screenVerticalMax,
                double dataHorizontalMin, double dataHorizontalMax, double dataVerticalMin, double dataVerticalMax) {
            this.screenHorizontalMin = screenHorizontalMin;
            this.screenHorizontalMax = screenHorizontalMax;
            this.screenVerticalMin = screenVerticalMin;
            this.screenVerticalMax = screenVerticalMax;

            this.dataHorizontalMin = dataHorizontalMin;
            this.dataHorizontalMax = dataHorizontalMax;
            this.dataVerticalMin = dataVerticalMin;
            this.dataVerticalMax = dataVerticalMax;

            this.horizontalResolution = (dataHorizontalMax - dataHorizontalMin) / (screenHorizontalMax - screenHorizontalMin);
            this.verticalResolution = (dataVerticalMax - dataVerticalMin) / (screenVerticalMax - screenVerticalMin);

        }
// note the vertical coordinates has been converted as we see usually

        /**
         *
         * @param Xa
         * @param Ya
         * @return
         */
        public double[] data2Screen(double Xa, double Ya) {
            double[] myout = new double[2];
            myout[0] = ((Xa - dataHorizontalMin) / horizontalResolution + screenHorizontalMin);
            myout[1] = (screenVerticalMax - (Ya - dataVerticalMin) / verticalResolution);
            return myout;
        }

        /**
         *
         * @param Xa
         * @param Ya
         * @return
         */
        public Point2D data2ScreenPoint(double Xa, double Ya) {
            Point2D myout = new Point2D.Double();
            myout.setLocation(((Xa - dataHorizontalMin) / horizontalResolution + screenHorizontalMin),
                    (screenVerticalMax - (Ya - dataVerticalMin) / verticalResolution));
            return myout;
        }

        /**
         *
         * @param X6
         * @param Y6
         * @return
         */
        public double[] screen2Data(int X6, int Y6) {
            double[] myout = new double[2];
            myout[0] = (X6 - screenHorizontalMin) * horizontalResolution + dataHorizontalMin;
            myout[1] = (screenVerticalMax - Y6) * verticalResolution + dataVerticalMin;
            return myout;
        }

        /**
         *
         * @param len
         * @return
         */
        public double horizontalSegmentData2Screen(double len) {
            return len / horizontalResolution;
        }

        /**
         *
         * @param len
         * @return
         */
        public double verticalSegmentData2Screen(double len) {
            return len / verticalResolution;
        }
    }

    /**
     *
     * @param xList
     * @param yList
     * @param title
     * @param outputPath
     * @throws Exception
     */
    public void drawScatterPlot(DoubleArrayList xList, DoubleArrayList yList, String title, String outputPath) throws Exception {
        BufferedImage image = new BufferedImage(this.canvasScope.width, this.canvasScope.height, BufferedImage.TYPE_INT_RGB);
        Graphics2D g2d = image.createGraphics();
        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        paintPlotArea(g2d);
        drawAxes(g2d, title, "p-value", "Probability");

        g2d.setStroke(new BasicStroke(1f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
        int pointSize = xList.size();
        //max and min value in vertical
        double vMax = Descriptive.max(yList);
        double vMin = Descriptive.min(yList);
        double hMin = 0;
        double hMax = 1;
        double thresholdForColor = 0.05;
        //sometimes there is only one point
        if (vMin - vMax == 0) {
            vMin = vMax - 1;
        }
        if (hMin - hMax == 0) {
            hMin = hMax - 1;
        }
        CoordinateTransformer transformer = new CoordinateTransformer();
        transformer.setupBasicScope(plotScope.x, plotScope.x + plotScope.width,
                plotScope.y, plotScope.y + plotScope.height, hMin, hMax, vMin, vMax);
        drawAxesScale(g2d, transformer);
        double rectangleSize = 4;
        for (int i = 0; i < pointSize; i++) {
            g2d.setColor(Color.BLACK);
            Point2D point1 = transformer.data2ScreenPoint(xList.getQuick(i), yList.getQuick(i));
            Rectangle2D.Double f = new Rectangle2D.Double(point1.getX() - rectangleSize / 2,
                    point1.getY() - rectangleSize / 2, rectangleSize, rectangleSize);
            g2d.draw(f);
            g2d.fill(f);
        }
        g2d.dispose();
        //outputJPEGFile(image, outputPath);
        outputPNGFile(image, outputPath);
    }

    /**
     *
     * @param xList
     * @param yList
     * @param title
     * @param outputPath
     * @throws Exception
     */
    public void drawHistogramPlot(DoubleArrayList xList, DoubleArrayList yList, String title, String outputPath) throws Exception {
        BufferedImage image = new BufferedImage(this.canvasScope.width, this.canvasScope.height, BufferedImage.TYPE_INT_RGB);
        Graphics2D g2d = image.createGraphics();
        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        paintPlotArea(g2d);
        drawAxes(g2d, title, "p-value", "Probability");

        g2d.setStroke(new BasicStroke(1f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
        int pointSize = xList.size();
        //max and min value in vertical
        double vMax = Descriptive.max(yList);
        double vMin = 0;
        double hMin = 0;
        double hMax = 1;
   
        //sometimes there is only one point
        if (vMin - vMax == 0) {
            vMin = vMax - 1;
        }
        if (hMin - hMax == 0) {
            hMin = hMax - 1;
        }
        CoordinateTransformer transformer = new CoordinateTransformer();
        transformer.setupBasicScope(plotScope.x, plotScope.x + plotScope.width,
                plotScope.y, plotScope.y + plotScope.height, hMin, hMax, vMin, vMax);
        drawAxesScale(g2d, transformer);
        Point2D corePoint = transformer.data2ScreenPoint(0, vMin);
        double rectangleSize = 1;
        g2d.setColor(Color.RED);
        for (int i = 0; i < pointSize; i++) {
            Point2D point1 = transformer.data2ScreenPoint(xList.getQuick(i), yList.getQuick(i));
            Rectangle2D.Double f = new Rectangle2D.Double(point1.getX() - rectangleSize / 2,
                    point1.getY(), rectangleSize, corePoint.getY() - point1.getY());
            g2d.draw(f);
            g2d.fill(f);
        }
        g2d.dispose();
        //outputJPEGFile(image, outputPath);
        outputPNGFile(image, outputPath);
    }
}
