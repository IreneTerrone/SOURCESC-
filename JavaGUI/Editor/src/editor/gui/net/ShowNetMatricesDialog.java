/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package editor.gui.net;

import editor.domain.NetObject;
import editor.domain.elements.GspnPage;
import editor.domain.unfolding.IncidenceMatrixFormatter;
import editor.domain.unfolding.MatrixMode;
import static editor.gui.net.NetEditorPanel.PAGE_BACKGROUND_COLOR;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.RenderingHints;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JComponent;
import javax.swing.Scrollable;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;
import latex.LatexFormula;

/**
 *
 * @author elvio
 */
public class ShowNetMatricesDialog extends javax.swing.JDialog {

    private final IncidenceMatrixFormatter matFormatter;
    LatexFormula latexFormula;
    LatexComponent latexComp;
    
    /**
     * Creates new form ShowNetMatricesDialog
     */
    public ShowNetMatricesDialog(java.awt.Frame parent, boolean modal, GspnPage gspn) {
        super(parent, modal);
        initComponents();
        
        latexComp = new LatexComponent();
        scrollPane.setViewportView(latexComp);
        scrollPane.setPreferredSize(new Dimension(600, 400));
                
        getRootPane().setDefaultButton(jButton_close);
        pack();
        setLocationRelativeTo(getOwner());
        
        matFormatter = new IncidenceMatrixFormatter(gspn);
        DefaultComboBoxModel<MatrixMode> model = new DefaultComboBoxModel<>();
        for (MatrixMode m : MatrixMode.values())
            model.addElement(m);
        model.setSelectedItem(MatrixMode.INCIDENCE_MATRIX);
        comboBox_matrixMode.setModel(model);
        
        update();
    }
    
    private void update() {
        System.out.println("update "+comboBox_matrixMode.getSelectedItem());
        String latex = matFormatter.latexFor((MatrixMode)comboBox_matrixMode.getSelectedItem(), 
                checkBox_showSumOfTerms.isSelected(), 
                checkBox_showZeros.isSelected(), 
                checkBox_showTermsStacked.isSelected());
        latexFormula = new LatexFormula(latex, NetObject.getUnitToPixels());
        
        latexComp.setPreferredSize(new Dimension(latexFormula.getWidth(), latexFormula.getHeight()));
        latexComp.repaint(); // drawing has changed
        latexComp.revalidate(); // area has changed
//        scrollPane.revalidate();
    }
    
    private class LatexComponent extends JComponent implements Scrollable {

        public LatexComponent() {
            setPreferredSize(new Dimension(400, 300));
        }
        

        @Override
        public void paint(Graphics g) {
            super.paint(g); //To change body of generated methods, choose Tools | Templates.
            System.out.println("paint!");
            
            Graphics2D g2 = (Graphics2D)g;
            g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
                                RenderingHints.VALUE_ANTIALIAS_ON);
            g2.setRenderingHint(RenderingHints.KEY_INTERPOLATION,
                                RenderingHints.VALUE_INTERPOLATION_BICUBIC);
            g2.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING,
                                RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
            g2.setRenderingHint(RenderingHints.KEY_FRACTIONALMETRICS,
                                RenderingHints.VALUE_FRACTIONALMETRICS_ON);

            // Draw the background
            g.setColor(PAGE_BACKGROUND_COLOR);
            g.fillRect(0, 0, getWidth(), getHeight());
            g.setColor(Color.BLACK);
            
            if (latexFormula != null) {
                latexFormula.draw(g2, 0, 0, 1.0, false);
            }
        }

        @Override
        public Dimension getPreferredScrollableViewportSize() {
            return getPreferredSize();
        }

        @Override
        public int getScrollableUnitIncrement(Rectangle visibleRect, int orientation, int direction) {
            return 12;
        }

        @Override
        public int getScrollableBlockIncrement(Rectangle visibleRect, int orientation, int direction) {
            if (orientation == SwingConstants.HORIZONTAL)
                return (int)visibleRect.getWidth();
            return (int)visibleRect.getHeight();
        }

        @Override
        public boolean getScrollableTracksViewportWidth() {
            return false;
        }

        @Override
        public boolean getScrollableTracksViewportHeight() {
            return false;
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
        java.awt.GridBagConstraints gridBagConstraints;

        southPanel = new javax.swing.JPanel();
        jButton_close = new javax.swing.JButton();
        jPanel1 = new javax.swing.JPanel();
        comboBox_matrixMode = new javax.swing.JComboBox<>();
        checkBox_showZeros = new javax.swing.JCheckBox();
        checkBox_showSumOfTerms = new javax.swing.JCheckBox();
        checkBox_showTermsStacked = new javax.swing.JCheckBox();
        scrollPane = new javax.swing.JScrollPane();

        setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);
        setTitle("Net matrices");

        southPanel.setLayout(new java.awt.GridBagLayout());

        jButton_close.setText("Close");
        jButton_close.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jButton_closeActionPerformed(evt);
            }
        });
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 0;
        gridBagConstraints.gridy = 1;
        gridBagConstraints.ipadx = 9;
        gridBagConstraints.insets = new java.awt.Insets(6, 6, 6, 6);
        southPanel.add(jButton_close, gridBagConstraints);

        jPanel1.setBorder(javax.swing.BorderFactory.createTitledBorder("Options"));
        jPanel1.setLayout(new java.awt.GridBagLayout());

        comboBox_matrixMode.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                comboBox_matrixModeActionPerformed(evt);
            }
        });
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.anchor = java.awt.GridBagConstraints.LINE_START;
        jPanel1.add(comboBox_matrixMode, gridBagConstraints);

        checkBox_showZeros.setSelected(true);
        checkBox_showZeros.setText("Show zeros.");
        checkBox_showZeros.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                checkBox_showZerosActionPerformed(evt);
            }
        });
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 0;
        gridBagConstraints.gridy = 1;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.LINE_START;
        gridBagConstraints.weightx = 1.0;
        jPanel1.add(checkBox_showZeros, gridBagConstraints);

        checkBox_showSumOfTerms.setSelected(true);
        checkBox_showSumOfTerms.setText("Show sum of terms.");
        checkBox_showSumOfTerms.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                checkBox_showSumOfTermsActionPerformed(evt);
            }
        });
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 0;
        gridBagConstraints.gridy = 2;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.LINE_START;
        gridBagConstraints.weightx = 1.0;
        jPanel1.add(checkBox_showSumOfTerms, gridBagConstraints);

        checkBox_showTermsStacked.setText("Show terms stacked");
        checkBox_showTermsStacked.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                checkBox_showTermsStackedActionPerformed(evt);
            }
        });
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 0;
        gridBagConstraints.gridy = 3;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.LINE_START;
        gridBagConstraints.weightx = 1.0;
        jPanel1.add(checkBox_showTermsStacked, gridBagConstraints);

        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.fill = java.awt.GridBagConstraints.BOTH;
        gridBagConstraints.weightx = 1.0;
        southPanel.add(jPanel1, gridBagConstraints);

        getContentPane().add(southPanel, java.awt.BorderLayout.SOUTH);
        getContentPane().add(scrollPane, java.awt.BorderLayout.CENTER);

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void jButton_closeActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jButton_closeActionPerformed
        SwingUtilities.invokeLater(new Runnable() {
            @Override
            public void run() {
                setVisible(false);
                dispose();
            }
        });
    }//GEN-LAST:event_jButton_closeActionPerformed

    private void comboBox_matrixModeActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_comboBox_matrixModeActionPerformed
        update();
    }//GEN-LAST:event_comboBox_matrixModeActionPerformed

    private void checkBox_showZerosActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_checkBox_showZerosActionPerformed
        update();
    }//GEN-LAST:event_checkBox_showZerosActionPerformed

    private void checkBox_showSumOfTermsActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_checkBox_showSumOfTermsActionPerformed
        update();
    }//GEN-LAST:event_checkBox_showSumOfTermsActionPerformed

    private void checkBox_showTermsStackedActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_checkBox_showTermsStackedActionPerformed
        update();
    }//GEN-LAST:event_checkBox_showTermsStackedActionPerformed

//    /**
//     * @param args the command line arguments
//     */
//    public static void main(String args[]) {
//        /* Set the Nimbus look and feel */
//        //<editor-fold defaultstate="collapsed" desc=" Look and feel setting code (optional) ">
//        /* If Nimbus (introduced in Java SE 6) is not available, stay with the default look and feel.
//         * For details see http://download.oracle.com/javase/tutorial/uiswing/lookandfeel/plaf.html 
//         */
//        try {
//            for (javax.swing.UIManager.LookAndFeelInfo info : javax.swing.UIManager.getInstalledLookAndFeels()) {
//                if ("Nimbus".equals(info.getName())) {
//                    javax.swing.UIManager.setLookAndFeel(info.getClassName());
//                    break;
//                }
//            }
//        } catch (ClassNotFoundException ex) {
//            java.util.logging.Logger.getLogger(ShowNetMatricesDialog.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
//        } catch (InstantiationException ex) {
//            java.util.logging.Logger.getLogger(ShowNetMatricesDialog.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
//        } catch (IllegalAccessException ex) {
//            java.util.logging.Logger.getLogger(ShowNetMatricesDialog.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
//        } catch (javax.swing.UnsupportedLookAndFeelException ex) {
//            java.util.logging.Logger.getLogger(ShowNetMatricesDialog.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
//        }
//        //</editor-fold>
//
//        /* Create and display the form */
//        java.awt.EventQueue.invokeLater(new Runnable() {
//            public void run() {
//                new ShowNetMatricesDialog().setVisible(true);
//            }
//        });
//    }

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JCheckBox checkBox_showSumOfTerms;
    private javax.swing.JCheckBox checkBox_showTermsStacked;
    private javax.swing.JCheckBox checkBox_showZeros;
    private javax.swing.JComboBox<MatrixMode> comboBox_matrixMode;
    private javax.swing.JButton jButton_close;
    private javax.swing.JPanel jPanel1;
    private javax.swing.JScrollPane scrollPane;
    private javax.swing.JPanel southPanel;
    // End of variables declaration//GEN-END:variables
}
