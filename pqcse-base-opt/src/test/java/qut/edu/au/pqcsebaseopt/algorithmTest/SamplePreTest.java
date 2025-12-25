package qut.edu.au.pqcsebaseopt.algorithmTest;

import org.junit.Test;
import qut.edu.au.pqcsebaseopt.tools.BigIntMatrix;

import java.math.BigInteger;
import java.security.SecureRandom;
import java.util.Arrays;

import static org.junit.Assert.assertArrayEquals;
import static qut.edu.au.pqcsebaseopt.algorithms.SamplePre.samplePre;
import static qut.edu.au.pqcsebaseopt.algorithms.SamplePre.solveLinearSystemFast;
import static qut.edu.au.pqcsebaseopt.algorithms.TrapGen.trapGen;

public class SamplePreTest {

    @Test
    public void testSamplePre() {
        BigInteger q = BigInteger.valueOf(17);
        int n = 2;
        int m = 4;

        // Define matrix A =
        // [1 2 3 4]
        // [0 1 5 6]
        BigIntMatrix A = new BigIntMatrix(n, m, q);
        A.set(0,0,1); A.set(0,1,2); A.set(0,2,3); A.set(0,3,4);
        A.set(1,0,0); A.set(1,1,1); A.set(1,2,5); A.set(1,3,6);

        // Target vector u = [5, 6]
        BigInteger[] u = {BigInteger.valueOf(5), BigInteger.valueOf(6)};

        System.out.println("Solving A*t = u...");
        BigInteger[] t = solveLinearSystemFast(A, u);

        System.out.println("Particular solution t: " + Arrays.toString(t));

        // Verify A*t = u
        boolean correct = true;
        for(int i=0; i<n; i++) {
            BigInteger sum = BigInteger.ZERO;
            for(int j=0; j<m; j++) {
                sum = sum.add(A.get(i, j).multiply(t[j]));
            }
            sum = sum.mod(q);
            if (!sum.equals(u[i])) correct = false;
        }
        System.out.println("Check A*t == u: " + correct);
    }

    @Test
    public void testSamplePre1() {
        double sigma = 100.0; // 高斯参数，需要足够大以平滑格
        BigInteger q = BigInteger.valueOf(17);
        int n = 5;
        BigIntMatrix[] ms = trapGen(n,q.longValue(),1);
        BigIntMatrix A = ms[0];
        BigIntMatrix T_A = ms[1];
        int m = A.getCols();
        BigInteger[] u = new BigInteger[n];
        SecureRandom rand = new SecureRandom();
        for(int i=0; i<n; i++) u[i] = new BigInteger(q.bitLength(), rand).mod(q);
        System.out.println("Target u: " + Arrays.toString(u));
        System.out.println("Solving A*t = u...");
        BigInteger[] t = solveLinearSystemFast(A, u);

        System.out.println("Particular solution t: " + Arrays.toString(t));

        BigInteger[] e = samplePre(A, T_A, sigma, u);

        System.out.println("Result e generated. Length: " + e.length);


        System.out.println("=== 4. Verification: A * e = u mod q ===");

        // 计算 A * e
        // 手动实现 Matrix * Vector
        BigInteger[] res = new BigInteger[n];
        for(int i=0; i<n; i++) {
            BigInteger sum = BigInteger.ZERO;
            for(int j=0; j<e.length; j++) {
                sum = sum.add(A.get(i, j).multiply(e[j]));
            }
            res[i] = sum.mod(q);
        }

        System.out.println("Calculated A*e: " + Arrays.toString(res));
        System.out.println("Expected u:     " + Arrays.toString(u));

        assertArrayEquals("SamplePre failed: A*e does not equal u!", u, res);

        System.out.println("Test Passed! SamplePre works correctly with TrapGen.");
    }
}
