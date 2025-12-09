package qut.edu.au.pqcsebase.algorithms;

import lombok.Getter;
import lombok.Setter;
import qut.edu.au.pqcsebase.tools.BigIntMatrix;
import qut.edu.au.pqcsebase.tools.HermiteNormalForm;
import java.security.SecureRandom;

import java.math.BigInteger;

public class TrapGenFirst {

    private static final SecureRandom secureRandom = new SecureRandom();

    private final long n_long;
    private final long q_long;
    // Use BigInteger parameters
    private final BigInteger n;
    private final BigInteger q;

    /**
     * m parameter (number of columns of A matrix) m = m1 + m2  A= [A1 | A2]
     */
    @Setter
    @Getter
    private int m;

    /**
     * m1 parameter (number of columns of A1 matrix)
     * m1 >= d = (1 + δ)n log2 q
     */
    @Setter@Getter
    private int m1;

    /**
     * m2 parameter (number of columns of A2 matrix)
     * m2 >= m1 * l  where l = ⎡log_r q⎤
     */
    @Setter@Getter
    private int m2;

    @Setter@Getter
    private int delta = 1; // δ > 0

    /**
     * r >= 2
     */
    @Setter@Getter
    private int r = 2;

    /**
     * l parameter (number of blocks) l = ⎡log_r q⎤
     */
    @Setter@Getter
    private  int l;

    public TrapGenFirst(long n, long q) {
        if (n <= 0) throw new IllegalArgumentException("n must be > 0");
        if (q <= 0) throw new IllegalArgumentException("q must be > 0");
        this.n_long = n;
        this.q_long = q;
        this.n = BigInteger.valueOf(n);
        this.q = BigInteger.valueOf(q);

        // Calculate parameters
        this.m1 = (1 + delta) * (int) n * (int) (Math.log(q) / Math.log(2));
        this.l = (int) Math.ceil(Math.log(q) / Math.log(r));
        this.m2 = l * m1;
        this.m = m1 + m2;
    }

    /**
     * Generate matrix A1 of size n x m1 with entries uniformly random in Z_q
     */
    public BigIntMatrix genA1() {
        BigIntMatrix A1 = BigIntMatrix.random((int)n_long, m1, q);
        // Simple rank check (Simplified: retry if image rows != n)
        if (BigIntMatrix.image(A1).getRowDimension() != n_long) {
            System.out.println("A1 not full rank, regenerating...");
            return genA1();
        }
        return A1;
    }

    /**
     * Generate HNF. Directly uses BigInteger[][] to prevent overflow.
     */
    public BigIntMatrix genH(BigIntMatrix matrix) {
        // Pass BigInteger[][] directly
        BigInteger[][] hnfRaw = HermiteNormalForm.hermiteNormalForm(matrix.getData());

        // Wrap back to BigIntMatrix
        BigIntMatrix H = new BigIntMatrix(hnfRaw.length, hnfRaw[0].length, q);
        for(int i=0; i<hnfRaw.length; i++) {
            for(int j=0; j<hnfRaw[0].length; j++) {
                H.set(i, j, hnfRaw[i][j]);
            }
        }
        return H;
    }

    /**
     * Generate matrix G.
     * Input matrix is H - I.
     * Fixes: Uses BigInteger division, no double casting, no modulo.
     */
    public BigIntMatrix genG(BigIntMatrix matrix) {
        BigIntMatrix G = new BigIntMatrix(m1, m2, q);
        //m1 block [G(1) | G(2) | ... | G(m1 - 1)| G(m1)] G(i) = [g_1(i), g_2(i), ..., g_l(i)]
        //G(i) m1 * l
        BigInteger bigR = BigInteger.valueOf(r);

        for (int i = 0; i < m1; i++) {
            for (int j = 0; j < l; j++) {
                // r_pow = r^(l - 1 - j)
                BigInteger rPow = bigR.pow(l - 1 - j);

                for (int k = 0; k < m1; k++) {
                    BigInteger val = matrix.get(k, i);

                    if (j != l - 1) {
                        // Correction: val / rPow (integer division), no modulo, no double
                        BigInteger scaled = val.divide(rPow);
                        G.set(k, i * l + j, scaled);
                    } else {
                        // Last column keeps original value (divide by 1)
                        G.set(k, i * l + j, val);
                    }
                }
            }
        }
        return G;
    }

    /**
     * Generate matrix P of size m2 x m1
     */
    public BigIntMatrix genP() {
        BigIntMatrix P = new BigIntMatrix(m2, m1, q);
        for (int j = 0; j < m1; j++) {
            int targetRowIndex = (j + 1) * l - 1;
            P.set(targetRowIndex, j, BigInteger.ONE);
        }
        return P;
    }

    /**
     * Generate Tl block.
     * Fixes: Removed double loop error. Only super-diagonal is -r.
     */
    public BigIntMatrix genTl(int l, BigInteger q) {
        BigIntMatrix Tl = BigIntMatrix.identity(l, q);
        BigInteger negR = BigInteger.valueOf(-r);

        // Only set super-diagonal (i, i+1)
        for (int i = 0; i < l - 1; i++) {
            Tl.set(i, i + 1, negR);
        }
        return Tl;
    }

    /**
     * Generate matrix U.
     */
    public BigIntMatrix genU() {
        // Initialize as Identity ensures bottom-right is I
        BigIntMatrix U = BigIntMatrix.identity(m2, q);
        BigIntMatrix Tl = genTl(l, q);

        for (int i = 0; i < m1; i++) {
            // Place Tl in diagonal blocks
            int startRow = i * l;
            int startCol = i * l;
            for (int r = 0; r < l; r++) {
                for (int c = 0; c < l; c++) {
                    U.set(startRow + r, startCol + c, Tl.get(r, c));
                }
            }
        }
        return U;
    }

    /**
     * Generate matrix R.
     */
    public BigIntMatrix genR() {
        BigIntMatrix R = new BigIntMatrix(m1, m2, q);
        double log2q = Math.log(q_long) / Math.log(2);
        long d = (long) Math.floor((1 + delta) * n_long * log2q);
        if (d > m1) d = m1;

        // Note: Java Random is not cryptographically secure, used for demo.
        for (int i = 0; i < d; i++) {
            for (int j = 0; j < m2; j++) {
                //rand = {0,1 |,2, |3}  0 (1/2), 1 (1/4), -1 (1/4)
                int rand = secureRandom.nextInt(4);
                if (rand <= 1) {
                    R.set(i, j, BigInteger.ZERO); // -1
                } else if (rand == 2) {
                    R.set(i, j, BigInteger.ONE); // 1
                } else {
                    R.set(i, j, BigInteger.valueOf(-1)); // 0
                }
            }
        }
        return R;
    }

    public void print(String name, BigIntMatrix matrix) {
        System.out.println(name + " dims: {" + matrix.getRowDimension() + ", " + matrix.getColumnDimension() + "}");
        // Avoid printing overly large matrices
        if (matrix.getRowDimension() <= 1000 && matrix.getColumnDimension() <= 1000) {
            System.out.println(name + ":\n" + matrix);
        } else {
            System.out.println(name + ": (Too large to print fully)");
        }
    }
}
