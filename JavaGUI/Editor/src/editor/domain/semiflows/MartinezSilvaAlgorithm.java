/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package editor.domain.semiflows;

import editor.domain.elements.Place;
import editor.domain.elements.Transition;
import editor.domain.struct.SilvaColom88Algorithm;
import java.util.ArrayList;
import java.util.Arrays;

/** Martinez-Silva version of the Farkas algorithm for the computation
 *  of the minimal set of P(T)-semiflows in a Petri net.
 *
 * @author elvio
 */
public class MartinezSilvaAlgorithm extends StructuralAlgorithm {

    
    // Iteration matrix [ D(i) | A(i) ], where D(0)=I and A(0) = Flow matrix
    // In addition, keep the B matrix for the Martinez-Silva optimization.
    // K starts with K=N. After the computation, K is the number of semiflows.
    public ArrayList<int[]> mD;     // KxN matrix
    public ArrayList<int[]> mA;     // KxM matrix
    public ArrayList<boolean[]> mB; // KxM matrix
    
    // This extra informations keep the initial place marking, when
    // the method is used for P-invariants. It allows to derive the place bounds
    // from the P-invariants
    public boolean computeBounds = false;
    public int[] initQuantity; // array of N elements
    public int[] lowerBnd, upperBnd; // place bounds
    
    
    SilvaColom88Algorithm scAlgo;
    

    // For P-semiflows: N=|P|, M=|T| (for T-semiflows: N=|T|, M=|P|)
    public MartinezSilvaAlgorithm(int N, int M) {
//        System.out.println("M="+M+", N="+N);
        super(N, M);
        mD = new ArrayList<>();
        mA = new ArrayList<>();
        mB = new ArrayList<>();
        for (int i = 0; i < N; i++) {
            mD.add(new int[N]);
            mA.add(new int[M]);
            mB.add(new boolean[M]);
            mD.get(i)[i] = 1;
        }
        scAlgo = new SilvaColom88Algorithm(N, M);
    }

    // Add a flow from i to j with the specified cardinality
    @Override
    public void addFlow(int i, int j, int card) {
//        System.out.println("addFlow("+i+", "+j+", "+card+");");
        mA.get(i)[j] += card;
        mB.get(i)[j] = (mA.get(i)[j] != 0);
        
        scAlgo.addFlow(i, j, card);
    }
    
    // Add the initial token of a place (only for P-invariant computation)
    @Override
    public void setInitQuantity(int i, int quantity) {
        assert i < N && i >= 0;
        if (initQuantity == null) { // Not initialized yet
            initQuantity = new int[N];
            computeBounds = true;
        }
        assert initQuantity[i] == 0;
        initQuantity[i] = quantity;
        
        scAlgo.setInitQuantity(i, quantity);
    }
    
    public int[] getSemiflow(int i) {
        return mD.get(i);
    }

    private static int gcd(int a, int b) {
        assert a >= 0 && b >= 0;
        if (a == 0)
            return b;

        while (b != 0) {
            if (a > b)
                a = a - b;
            else
                b = b - a;
        }

        return a;
    }
    
    private static int sign(int num) {
        if (num > 0)
            return +1;
        else if (num < 0)
            return -1;
        return 0;
    }
    

    @Override
    public void compute(boolean log, ProgressObserver obs) throws InterruptedException {
        if (log)
            System.out.println(this);
        // Matrix A starts with the flow matrix, D is the identity.
        // for every transition i=[0,M), repeat:
        for (int i = 0; i < M; i++) {
            if (log)
                System.out.println("\nStep "+i+"/"+(M-1)+"\n"+this);
            // Append to the matrix [D|A] every rows resulting as a non-negative
            // linear combination of row pairs from [D|A] whose sum zeroes
            // the i-th column of A.
            int nRows = numSemiflows();
            for (int r1 = 0; r1 < nRows; r1++) {
                if (mA.get(r1)[i] == 0) {
                    continue;
                }
                obs.advance(i, M+1, r1, nRows);
                for (int r2 = r1 + 1; r2 < nRows; r2++) {
                    checkInterrupted();
                    // Find two rows r1 and r2 such that r1[i] and r2[i] have opposite signs.
                    if (mA.get(r2)[i] == 0)
                        continue;
                    if (sign(mA.get(r1)[i]) == sign(mA.get(r2)[i]))
                        continue;
                    int abs1 = Math.abs(mA.get(r1)[i]), abs2 = Math.abs(mA.get(r2)[i]);

                    // Create a new row nr' such that:
                    //   nr = |r2[i]| * r1 + |ri[i]| * r2
                    //   nr' = nr / gcd(nr)
                    //   nr(B) = logical union of B[r1] and B[r2]
                    int[] nrA = new int[M];
                    int[] nrD = new int[N];
                    boolean[] nrB = new boolean[M];
                    int gcdAD = -1;
                    for (int k = 0; k < M; k++) {
                        nrA[k] = abs2 * mA.get(r1)[k] + abs1 * mA.get(r2)[k];
                        gcdAD = (k == 0) ? Math.abs(nrA[k]) : gcd(gcdAD, Math.abs(nrA[k]));
                        nrB[k] = mB.get(r1)[k] || mB.get(r2)[k];
                    }
                    assert nrA[i] == 0;
                    for (int k = 0; k < N; k++) {
                        nrD[k] = abs2 * mD.get(r1)[k] + abs1 * mD.get(r2)[k];
                        gcdAD = gcd(gcdAD, Math.abs(nrD[k]));
                    }
                    if (gcdAD != 1) {
//                        System.out.println("  gcdAD = " + gcdAD);
                        for (int k = 0; k < M; k++) {
                            nrA[k] /= gcdAD;
                        }
                        for (int k = 0; k < N; k++) {
                            nrD[k] /= gcdAD;
                        }
                    }
                    int nnzD = 0, ntrueB = 0;
                    for (int k = 0; k < M; k++) {
                        ntrueB += (nrB[k] ? 1 : 0);
                    }
                    for (int k = 0; k < N; k++) {
                        nnzD += (nrD[k] != 0 ? 1 : 0);
                    }

                    // Martinez-Silva optimization of the Farkas algorithm.
                    // The row is not a minimal support if the count of non-zero entries in D
                    // is greater than the number of TRUEs in B (for that row) + 1.
                    // If this happen, the row is not a minimal P(T)-semiflow and
                    // can be safely discarded.
                    if (nnzD > ntrueB+1) {
                        continue;
                    }
                    if (log)
                        System.out.println(i + ": ADD row " + r1 + " + row " + r2 +
                                           "  nnz(D)=" + nnzD + " r'=" + ntrueB);

                    mA.add(nrA);
                    mD.add(nrD);
                    mB.add(nrB);
                }
            }
            checkInterrupted();

            // Eliminate from [D|A] the rows in which the i-th column of A is not zero.
            int rr = nRows;
            while (rr > 0) {
                rr--;
                obs.advance(i, M+1, rr, nRows);
                if (mA.get(rr)[i] == 0) {
                    continue;
                }
                if (log)
                    System.out.println(i + ": DEL row " + rr);
                mA.remove(rr);
                mD.remove(rr);
                mB.remove(rr);
            }
            
            // Eliminate frm [D|A] the rows that are not minimal, doing an exhaustive search.
            removeNonMinimalSemiflows(log, obs);
        }
        if (log)
            System.out.println("\nRESULT:\n"+this);

        obs.advance(M+1, M+1, 1, 1);
        
        if (computeBounds) {
            computeBoundsFromInvariants();
            scAlgo.compute(log, obs);
        }
        
        setComputed();
    }
    
    // Test all semiflows exaustively to check if they are minimal
    // To do so, we take every pair of semiflows and we test if one is a
    // linear component of another. If it is so, subtract the component semiflow
    // from the other. When a semiflow becomes zero, it is removed from A,D and B.
    private void removeNonMinimalSemiflows(boolean log, ProgressObserver obs) 
            throws InterruptedException 
    {
        int rr = numSemiflows();
        while (rr > 0) {
            rr--;
            obs.advance(M, M+1, numSemiflows()-rr, numSemiflows());
            for (int i=0; i<numSemiflows(); i++) {
                checkInterrupted();
                if (i == rr)
                    continue;
                // Check if the semiflow D[rr] contains D[i]
                int mult = -1;
                boolean isComp = true; // Says if i is a linear component of rr
                for (int k=0; k<N; k++) {
                    if (mD.get(i)[k] != 0) {
                        if (mD.get(i)[k] > mD.get(rr)[k]) {
                            isComp = false;
                            break; // i is not a linear component of rr
                        }
                        if (mult == -1)
                            mult = mD.get(rr)[k] / mD.get(i)[k];
                        else
                            mult = Math.min(mult, mD.get(rr)[k] / mD.get(i)[k]);
                    }
                }
                if (!isComp)
                    continue;
                // in matrix D:  D[rr] -= mult * D[i]
                assert mult > 0;
                boolean isZero = true;
                for (int k=0; k<N; k++) {
                    mD.get(rr)[k] -= mult * mD.get(i)[k];
                    isZero = isZero && (mD.get(rr)[k] == 0);
                }
                if (isZero) {
                    // rr was a linear combination of other semiflows. Remove it
                    if (log)
                        System.out.println("DEL row " + rr);
                    mA.remove(rr);
                    mD.remove(rr);
                    mB.remove(rr);
                    break;
                }
            }
        }        
    }

    public int numSemiflows() {
        return mD.size();
    }
    
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < numSemiflows(); i++) {
            sb.append(i < 10 ? " " : "").append(i).append(":  ");
            for (int j = 0; j < N; j++) {
                sb.append(mD.get(i)[j] < 0 ? "" : " ").append(mD.get(i)[j]).append(" ");
            }
            sb.append("| ");
            for (int j = 0; j < M; j++) {
                sb.append(mA.get(i)[j] < 0 ? "" : " ").append(mA.get(i)[j]).append(" ");
            }
            sb.append("|");
            for (int j = 0; j < M; j++) {
                sb.append(mB.get(i)[j] ? " T" : " .");
            }
            sb.append("\n");
        }
        return sb.toString();
    }
    
    
    void computeBoundsFromInvariants() {
        lowerBnd = new int[N];
        upperBnd = new int[N];
        Arrays.fill(upperBnd, Integer.MAX_VALUE);
        
        // Read all place invariants
        for (int i=0; i<numSemiflows(); i++) {
            int[] inv = getSemiflow(i);
            int tokenCnt = 0;
            for (int p=0; p<N; p++) {
                tokenCnt += inv[p] * initQuantity[p];
            }
            
//            System.out.print("PINV: ");
//            for (int p=0; p<N; p++)
//                if (inv[p] > 0)
//                    System.out.print(inv[p]+"*P"+p+"(m0="+initQuantity[p]+")  ");
//            System.out.println("   tc="+tokenCnt);
            
            int num_nnz = 0, kk = -1, last_p = -1;
            for (int p=0; p<N; p++) {
                if (inv[p] > 0) {
                    last_p = p;
                    num_nnz++;
                    kk = tokenCnt / inv[p];
                    upperBnd[p] = Math.min(upperBnd[p], kk);
                }
            }
            if (num_nnz == 1 && kk > lowerBnd[last_p])
                lowerBnd[last_p] = kk;
        }
        
//        for (int p=0; p<N; p++) {
//            System.out.println("Place "+p+" has bounds: ["+lowerBnd[p]+" - "+upperBnd[p]+"]");
//        }
    }
    
    public int getUpperBoundOf(int p) {
        assert isComputed() && computeBounds && initQuantity != null && p >= 0 && p < N;
        return upperBnd[p];
    }
    
    public int getLowerBoundOf(int p) {
        assert isComputed() && computeBounds && initQuantity != null && p >= 0 && p < N;
        return lowerBnd[p];
    }
    
    
    
    
    public String toLatexString(SemiFlows.Type type, ArrayList<Place> places, 
                                ArrayList<Transition> transitions,
                                boolean showZeros) 
    {
        boolean add_m0_col = (type==SemiFlows.Type.PLACE_SEMIFLOW);
        StringBuilder sb = new StringBuilder();
        
        // header
        sb.append("$\\begin{array}{r");
        for (int f=0; f<numSemiflows(); f++)
            sb.append("|c");
        if (add_m0_col)
            sb.append("|r");
        sb.append("}\n ");
        for (int f=0; f<numSemiflows(); f++)
            sb.append("& i_{").append(f+1).append("}");
        if (add_m0_col)
            sb.append("& \\mathbf{m}_0");
        sb.append("\\\\ \n\\hline\n");
                
        // row for a place/transition
        for (int pl=0; pl<N; pl++) {
//            String m0p = null;
//            if (add_m0_col) {
//                m0p = places.get(pl).getInitMarkingExpr();
//                if (m0p.isBlank())
//                    m0p = "0";
//                int m0val = -1000;
//                try {
//                    m0val = Integer.parseInt(m0p);
//                }
//                catch (NumberFormatException e) {}
//            }
            
            if (type == SemiFlows.Type.PLACE_SEMIFLOW)
                sb.append(places.get(pl).getUniqueNameDecor().getLatexFormula().getLatex());
            else
                sb.append(transitions.get(pl).getUniqueNameDecor().getLatexFormula().getLatex());
            
            for (int f=0; f<numSemiflows(); f++) {
                int[] semiflow = getSemiflow(f);
                
                sb.append(" &");
                if (showZeros || semiflow[pl]!=0) {
                    String color;
                    if (semiflow[pl] > 0)
                        color = "Blue";
                    else if (semiflow[pl] < 0)
                        color = "Mahogany";
                    else
                        color = "Gray";
                    sb.append("\\textcolor{").append(color).append("}{").append(semiflow[pl]).append("}");
                }
            } 
            
            if (add_m0_col) {
                sb.append(" & \\mathbf{").append(initQuantity[pl]).append("}");
            }
            sb.append("\\\\ \n\\hline\n");
            
        }
        
        // final row
        if (type == SemiFlows.Type.PLACE_SEMIFLOW) {
            sb.append("\n \\mathbf{m}_0 \\cdot I & ");
            
            for (int f=0; f<numSemiflows(); f++) {
                int[] semiflow = getSemiflow(f);
                int sum = 0;
                for (int j=0; j<semiflow.length; j++) {
                    int initMark = initQuantity[j];
                    sum += initMark * semiflow[j];
                }
                
                sb.append("\\mathbf{").append(sum).append("} & ");
            }
            sb.append("\\\\ \n");
        }
//        sb.append("\n& ");
//        for (int f=0; f<pinMat.size(); f++) {
//            if (!extra_pm0[f].isEmpty()) {
//                if (sum_pm0[f] == 0)
//                    sb.append(extra_pm0[f]);
//                else
//                    sb.append(sum_pm0[f]+"+"+extra_pm0[f]);
//            }
//            else sb.append(sum_pm0[f]);
//            sb.append("& ");
//        }
//        sb.append("\\\\ \n");
        
        
        sb.append("\\end{array}$");
        return sb.toString();
    }

//    public static void main(String[] args) {
//        int NP = 14, MT = 10;
//        MartinezSilvaAlgorithm msa = new MartinezSilvaAlgorithm(NP, MT);
//        int[][] flow = {
//            {5, 4}, {1, 3}, // t1
//            {14, 7}, {6, 8}, // t2
//            {9}, {10, 11}, // t3
//            {2, 13}, {4}, // t4
//            {1}, {2}, // t5
//            {3}, {14}, // t6
//            {6}, {5}, // t7
//            {8}, {9}, // t8
//            {10, 12}, {7}, // t9
//            {11}, {12, 13} // t10
//        };
//        for (int t = 0; t < MT; t++) {
//            int[] in = flow[2 * t], out = flow[2 * t + 1];
//            for (int p = 0; p < in.length; p++) {
//                msa.addFlow(in[p] - 1, t, -1);
//            }
//            for (int p = 0; p < out.length; p++) {
//                msa.addFlow(out[p] - 1, t, +1);
//            }
//        }
//        msa.compute(true);
//    }
    
    public static MartinezSilvaAlgorithm init1() {
        int M=16, N=12;
        MartinezSilvaAlgorithm msa = new MartinezSilvaAlgorithm(N, M);
        msa.addFlow(1, 1, 1);
        msa.addFlow(0, 2, 1);
        msa.addFlow(2, 3, 1);
        msa.addFlow(2, 1, 1);
        msa.addFlow(1, 0, 1);
        msa.addFlow(0, 0, 1);
        msa.addFlow(3, 2, 1);
        msa.addFlow(3, 3, 1);
        msa.addFlow(1, 5, 1);
        msa.addFlow(8, 5, -1);
        msa.addFlow(8, 6, 1);
        msa.addFlow(4, 6, -1);
        msa.addFlow(4, 4, 1);
        msa.addFlow(1, 4, -1);
        msa.addFlow(4, 1, -1);
        msa.addFlow(7, 1, -1);
        msa.addFlow(5, 2, -1);
        msa.addFlow(6, 2, -1);
        msa.addFlow(0, 7, -1);
        msa.addFlow(0, 8, 1);
        msa.addFlow(6, 8, -1);
        msa.addFlow(6, 9, 1);
        msa.addFlow(9, 9, -1);
        msa.addFlow(9, 7, 1);
        msa.addFlow(3, 13, 1);
        msa.addFlow(10, 13, -1);
        msa.addFlow(10, 14, 1);
        msa.addFlow(5, 14, -1);
        msa.addFlow(5, 10, 1);
        msa.addFlow(3, 10, -1);
        msa.addFlow(2, 12, -1);
        msa.addFlow(2, 11, 1);
        msa.addFlow(7, 11, -1);
        msa.addFlow(7, 15, 1);
        msa.addFlow(11, 15, -1);
        msa.addFlow(11, 12, 1);
        msa.addFlow(8, 0, -1);
        msa.addFlow(9, 0, -1);
        msa.addFlow(11, 3, -1);
        msa.addFlow(10, 3, -1);
        return msa;
    }
    

    public static void main(String[] args) throws InterruptedException {
        MartinezSilvaAlgorithm msa = init1();
        ProgressObserver obs = new ProgressObserver() {
            @Override
            public void advance(int step, int total, int s, int t) {
            }
        };
        msa.compute(true, obs);
    }
}
