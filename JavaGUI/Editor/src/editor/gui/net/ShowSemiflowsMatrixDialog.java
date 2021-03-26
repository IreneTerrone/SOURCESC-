/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package editor.gui.net;

import editor.domain.NetObject;
import editor.domain.semiflows.MartinezSilvaAlgorithm;
import editor.domain.semiflows.NetIndex;
import editor.domain.semiflows.SemiFlows;
import editor.gui.ResourceFactory;
import java.awt.Dimension;
import javax.swing.SwingUtilities;
import latex.JLatexComponent;
import latex.LatexFormula;

/**
 *
 * @author elvio
 */
public class ShowSemiflowsMatrixDialog extends javax.swing.JDialog {

    private final SemiFlows.Type type;
    private final String netName;
    private final MartinezSilvaAlgorithm algo;
    private final NetIndex netIndex;
    JLatexComponent latexComp;
    
    /**
     * Creates new form ShowNetMatricesDialog
     */
    public ShowSemiflowsMatrixDialog(java.awt.Frame parent, boolean modal, MartinezSilvaAlgorithm algo,
                                     SemiFlows.Type type, String netName, NetIndex netIndex) {
        super(parent, modal);
        this.algo = algo;
        this.type = type;
        this.netName = netName;
        this.netIndex = netIndex;
        initComponents();
        
        latexComp = new JLatexComponent();
        scrollPane.setViewportView(latexComp);
        scrollPane.setPreferredSize(new Dimension(600, 400));
        
        button_copyToClipboard.setIcon(ResourceFactory.getInstance().getPageWhitePaste16());
        
        getRootPane().setDefaultButton(jButton_close);
        pack();
        setLocationRelativeTo(getOwner());
        
        update();
    }
    
    private void update() {
        String latex = algo.toLatexString(type, netIndex, checkBox_showZeros.isSelected());
        LatexFormula formula = new LatexFormula(latex, NetObject.getUnitToPixels());
        
        latexComp.setFormula(formula);
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
        checkBox_showZeros = new javax.swing.JCheckBox();
        button_copyToClipboard = new javax.swing.JButton();
        button_saveAsPdf = new javax.swing.JButton();
        scrollPane = new javax.swing.JScrollPane();

        setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);
        setTitle("Semiflows matrix");

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

        button_copyToClipboard.setText("CopyToClipboard");
        button_copyToClipboard.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                button_copyToClipboardActionPerformed(evt);
            }
        });
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 1;
        gridBagConstraints.gridy = 3;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.LINE_END;
        jPanel1.add(button_copyToClipboard, gridBagConstraints);

        button_saveAsPdf.setText("Save as PDF...");
        button_saveAsPdf.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                button_saveAsPdfActionPerformed(evt);
            }
        });
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 1;
        gridBagConstraints.gridy = 2;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        jPanel1.add(button_saveAsPdf, gridBagConstraints);

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

    private void checkBox_showZerosActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_checkBox_showZerosActionPerformed
        update();
    }//GEN-LAST:event_checkBox_showZerosActionPerformed

    private void button_copyToClipboardActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_button_copyToClipboardActionPerformed
        latexComp.copyToClipboard();
    }//GEN-LAST:event_button_copyToClipboardActionPerformed

    private void button_saveAsPdfActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_button_saveAsPdfActionPerformed
        latexComp.saveAsPdf(type.printableName()+" of "+netName);
    }//GEN-LAST:event_button_saveAsPdfActionPerformed



    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton button_copyToClipboard;
    private javax.swing.JButton button_saveAsPdf;
    private javax.swing.JCheckBox checkBox_showZeros;
    private javax.swing.JButton jButton_close;
    private javax.swing.JPanel jPanel1;
    private javax.swing.JScrollPane scrollPane;
    private javax.swing.JPanel southPanel;
    // End of variables declaration//GEN-END:variables
}
