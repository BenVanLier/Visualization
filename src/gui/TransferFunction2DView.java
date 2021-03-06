/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package gui;

import java.awt.Color;
import java.awt.Cursor;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionAdapter;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Rectangle2D;

/**
 *
 * @author michel
 */
public class TransferFunction2DView extends javax.swing.JPanel {

    TransferFunction2DEditor ed;
    private final int DOTSIZE = 8;
    public Ellipse2D.Double baseControlPoint, radiusControlPoint;
    boolean selectedBaseControlPoint, selectedRadiusControlPoint;
    private double maxGradientMagnitude;
    
    /**
     * Creates new form TransferFunction2DView
     * @param ed
     */
    public TransferFunction2DView(TransferFunction2DEditor ed) {
        initComponents();
        
        this.ed = ed;
        maxGradientMagnitude = ed.getMGM();
        selectedBaseControlPoint = false;
        selectedRadiusControlPoint = false;
        addMouseMotionListener(new TriangleWidgetHandler());
        addMouseListener(new SelectionHandler());
    }
    
    @Override
    public void paintComponent(Graphics g) {

        Graphics2D g2 = (Graphics2D) g;

        int w = this.getWidth();
        int h = this.getHeight();
        g2.setColor(Color.white);
        g2.fillRect(0, 0, w, h);
        
        double maxHistoMagnitude = ed.histogram[0];
        for (int i = 0; i < ed.histogram.length; i++) {
            maxHistoMagnitude = ed.histogram[i] > maxHistoMagnitude ? ed.histogram[i] : maxHistoMagnitude;
        }
        
        double binWidth = (double) w / (double) ed.xbins;
        double binHeight = (double) h / (double) ed.ybins;
        maxHistoMagnitude = Math.log(maxHistoMagnitude);
        
        
        for (int y = 0; y < ed.ybins; y++) {
            for (int x = 0; x < ed.xbins; x++) {
                if (ed.histogram[y * ed.xbins + x] > 0) {
                    int intensity = (int) Math.floor(255 * (1.0 - Math.log(ed.histogram[y * ed.xbins + x]) / maxHistoMagnitude));
                    g2.setColor(new Color(intensity, intensity, intensity));
                    g2.fill(new Rectangle2D.Double(x * binWidth, h - (y * binHeight), binWidth, binHeight));
                }
            }
        }

        int yposB = (int) Math.floor((ed.ybins / maxGradientMagnitude) * ed.triangleWidget.minMagDraw);
        int yposR = (int) Math.floor((ed.ybins / maxGradientMagnitude) * ed.triangleWidget.maxMagDraw);
        int xpos = (int) (ed.triangleWidget.baseIntensity * binWidth);
        
        g2.setColor(Color.black);
        baseControlPoint = new Ellipse2D.Double(xpos - DOTSIZE / 2, yposB - DOTSIZE, DOTSIZE, DOTSIZE);
        g2.fill(baseControlPoint);
        g2.drawLine(xpos, yposB, xpos - (int) (ed.triangleWidget.radius * binWidth * maxGradientMagnitude), yposR);
        g2.drawLine(xpos, yposB, xpos + (int) (ed.triangleWidget.radius * binWidth * maxGradientMagnitude), yposR);
        g2.drawLine(xpos - (int) (ed.triangleWidget.radius * binWidth * maxGradientMagnitude), yposR, xpos + (int) (ed.triangleWidget.radius * binWidth * maxGradientMagnitude), yposR);
        radiusControlPoint = new Ellipse2D.Double(xpos + (ed.triangleWidget.radius * binWidth * maxGradientMagnitude) - DOTSIZE / 2,  yposR, DOTSIZE, DOTSIZE);
        g2.fill(radiusControlPoint);
    }
    
    
    private class TriangleWidgetHandler extends MouseMotionAdapter {

        @Override
        public void mouseMoved(MouseEvent e) {
            if (baseControlPoint.contains(e.getPoint()) || radiusControlPoint.contains(e.getPoint())) {
                setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
            } else {
                setCursor(Cursor.getDefaultCursor());
            }
        }
        
        @Override
        public void mouseDragged(MouseEvent e) {
            if (selectedBaseControlPoint || selectedRadiusControlPoint) {
                Point dragEnd = e.getPoint();
                
                if (selectedBaseControlPoint) {
                    // restrain to horizontal movement
                    // also make sure that baseControlPoint cannot go above the radiusControlpoint
                    dragEnd.setLocation(dragEnd.x, dragEnd.y);
                    if (dragEnd.y < radiusControlPoint.getCenterY()) {
                        dragEnd.y = (int) radiusControlPoint.getCenterY();
                    }
                } else if (selectedRadiusControlPoint) {
                    // restrain to horizontal movement and avoid radius getting 0
                    // also make sure that radiusControlPoint cannot go under the baseControlpoint
                    dragEnd.setLocation(dragEnd.x, dragEnd.y);
                    if (dragEnd.x - baseControlPoint.getCenterX() <= 0) {
                        dragEnd.x = (int) (baseControlPoint.getCenterX() + 1);
                    }
                    if (dragEnd.y > baseControlPoint.getCenterY()) {
                        dragEnd.y = (int) (baseControlPoint.getCenterY());
                    }
                }
                if (dragEnd.x < 0) {
                    dragEnd.x = 0;
                }
                if (dragEnd.x >= getWidth()) {
                    dragEnd.x = getWidth() - 1;
                }
                if (dragEnd.y < 0) {
                    dragEnd.y = 0;
                }
                if (dragEnd.y >= getHeight()) {
                    dragEnd.y = getHeight() - 1;
                }
                double w = getWidth();
                double h = getHeight();
                double binWidth = (double) w / (double) ed.xbins;
                double binHeight = (double) h / (double) ed.ybins;
                if (selectedBaseControlPoint) {
                    ed.triangleWidget.baseIntensity = (short) (dragEnd.x / binWidth);
                    ed.triangleWidget.minMagDraw = (dragEnd.y * binHeight *(maxGradientMagnitude/ed.ybins));
                    ed.triangleWidget.minMag = maxGradientMagnitude - ed.triangleWidget.minMagDraw;
                } else if (selectedRadiusControlPoint) {
                    ed.triangleWidget.radius = (dragEnd.x - (ed.triangleWidget.baseIntensity * binWidth))/(binWidth * maxGradientMagnitude);
                    ed.triangleWidget.maxMagDraw = (dragEnd.y * binHeight *(maxGradientMagnitude/ed.ybins));
                    ed.triangleWidget.maxMag = maxGradientMagnitude - ed.triangleWidget.maxMagDraw;
                }
                ed.setSelectedInfo();
                
                repaint();
            } 
        }

    }
    
    
    private class SelectionHandler extends MouseAdapter {
        @Override
        public void mousePressed(MouseEvent e) {
            if (baseControlPoint.contains(e.getPoint())) {
                selectedBaseControlPoint = true;
            } else if (radiusControlPoint.contains(e.getPoint())) {
                selectedRadiusControlPoint = true;
            } else {
                selectedRadiusControlPoint = false;
                selectedBaseControlPoint = false;
            }
        }
        
        @Override
        public void mouseReleased(MouseEvent e) {
            selectedRadiusControlPoint = false;
            selectedBaseControlPoint = false;
            ed.changed();
            repaint();
        }
    }
    
    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(this);
        this.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 400, Short.MAX_VALUE)
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 300, Short.MAX_VALUE)
        );
    }// </editor-fold>//GEN-END:initComponents


    // Variables declaration - do not modify//GEN-BEGIN:variables
    // End of variables declaration//GEN-END:variables
}
