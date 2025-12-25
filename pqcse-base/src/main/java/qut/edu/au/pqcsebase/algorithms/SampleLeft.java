package qut.edu.au.pqcsebase.algorithms;

import qut.edu.au.pqcsebase.tools.BigIntMatrix;

import java.math.BigInteger;
import java.util.Arrays;

import static qut.edu.au.pqcsebase.algorithms.SamplePre.samplePre;
import static qut.edu.au.pqcsebase.algorithms.SampleZ.sampleZ;

public class SampleLeft {

    /**
     * @param A  The left matrix (n x m1) which has a trapdoor.  [A|M1]e = u  n*(m1+m2) (m1+m2)*1 n*1
     * @param M1 The right matrix (n x m2) which is arbitrary/random.
     * @param TA the trapdoor of A
     * @param sigma  Gaussian parameter
     * @param u  target vector (length n)
     * @param nParam security parameter
     * @return e e = [e1 ; e3] A short vector e = (e1, e2) such that A*e1 + M*e2 = u mod q.
     */
    public static BigInteger[] sampleLeft(BigIntMatrix A, BigIntMatrix M1, BigIntMatrix TA, double sigma, BigInteger[] u, int nParam) {

        int n = A.getRowDimension();
        int m_left = A.getColumnDimension();
        int m_right = M1.getColumnDimension();
        BigInteger q = A.getModulus();

        //1. sample the random part e2
        BigInteger[] e2 = new BigInteger[m_right];
        for (int i = 0; i < m_right; i++) {
            e2[i] = sampleZ(sigma, 0.0, nParam);
        }

        //2. Compute the shifted target y = A*e1  v = M*e2 n*1
        // solve: A*e1 + M*e2 = u  ->  A*e1 = u - M1*e2 mod q   --> y = u - v
        BigInteger[] v = new BigInteger[n];
        Arrays.fill(v, BigInteger.ZERO);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m_right; j++) {
                // v[i] += M[i][j] * e2[j]
                BigInteger val = M1.get(i, j).multiply(e2[j]);
                v[i] = v[i].add(val);
            }
            v[i] = v[i].mod(q);
        }

        BigInteger[] y = new BigInteger[n];
        for (int i = 0; i < n; i++) {
            y[i] = u[i].subtract(v[i]).mod(q);
        }

        //3.Sample the trapdoor part e1  A*e1 = y  compute e1 use samplePre
        BigInteger[] e1 = samplePre(A,TA,sigma,y);

        BigInteger[] e = new BigInteger[m_left + m_right];
        System.arraycopy(e1, 0, e, 0, m_left);
        System.arraycopy(e2, 0, e, m_left, m_right);

        return e;
    }
}
