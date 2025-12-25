package qut.edu.au.pqcsebaseopt.algorithms;


import qut.edu.au.pqcsebaseopt.tools.BigIntMatrix;

import java.math.BigInteger;
import java.util.stream.IntStream;

import static qut.edu.au.pqcsebaseopt.algorithms.SamplePre.samplePre;
import static qut.edu.au.pqcsebaseopt.algorithms.SampleZ.sampleZ;

public class SampleLeft {

    public static BigInteger[] sampleLeft(BigIntMatrix A, BigIntMatrix M1, BigIntMatrix TA, double sigma, BigInteger[] u, int nParam) {

        int n = A.getRows();
        int m_left = A.getCols();
        int m_right = M1.getCols();
        BigInteger q = A.getModulus();

        // 1. Sample e2 (Parallel sampling)
        BigInteger[] e2 = new BigInteger[m_right];
        IntStream.range(0, m_right).parallel().forEach(i -> {
            e2[i] = sampleZ(sigma, 0.0, nParam);
        });

        // 2. Compute v = M1 * e2 (Optimized Matrix-Vector Multiplication)
        // Uses the optimized parallel implementation in BigIntMatrix
        BigInteger[] v = BigIntMatrix.multiplyMatrixVector(M1, e2);

        // Compute y = u - v mod q
        BigInteger[] y = new BigInteger[n];
        IntStream.range(0, n).parallel().forEach(i -> {
            y[i] = u[i].subtract(v[i]).mod(q);
        });

        // 3. Sample e1 using SamplePre
        BigInteger[] e1 = samplePre(A, TA, sigma, y);

        // 4. Combine result
        BigInteger[] e = new BigInteger[m_left + m_right];
        System.arraycopy(e1, 0, e, 0, m_left);
        System.arraycopy(e2, 0, e, m_left, m_right);

        return e;
    }
}
