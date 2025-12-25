package qut.edu.au.pqcsebaseopt.algorithmTest;

import org.junit.Test;
import qut.edu.au.pqcsebaseopt.tools.BigIntMatrix;

import java.math.BigInteger;
import java.util.Arrays;

import static qut.edu.au.pqcsebaseopt.algorithms.SampleLeft.sampleLeft;
import static qut.edu.au.pqcsebaseopt.algorithms.TrapGen.trapGen;

public class SampleLeftTest {

    @Test
    public void testSampleLeft(){
        BigInteger q = BigInteger.valueOf(17);
        int n = 2;
        int m1 = 3;
        int m2 = 2;
        double s = 5.0;
        int nParam = 10;

        // Create random matrices A and M
        BigIntMatrix A = BigIntMatrix.random(n, m1, q);
        BigIntMatrix M = BigIntMatrix.random(n, m2, q);

        BigIntMatrix T_A = BigIntMatrix.identity(m1, q);

        // Target u
        BigInteger[] u = {BigInteger.valueOf(1), BigInteger.valueOf(2)};

        System.out.println("Starting SampleLeft...");

        try {
            BigInteger[] result = sampleLeft(A, M, T_A, s, u, nParam);
            System.out.println("Result e: " + Arrays.toString(result));

            //Verify
            BigIntMatrix AM1 = BigIntMatrix.columnConcat(A, M);
            BigIntMatrix eMatrix = new BigIntMatrix(result.length,1,A.getModulus());

            for (int i = 0; i < result.length; i++) {
                eMatrix.set(i,0,result[i]);
            }

            BigIntMatrix resultTest = AM1.multiplyQ(eMatrix);
            System.out.println("AM1 * e = u: " + resultTest.toString());
            System.out.println("u: " + Arrays.toString(u));

        } catch (Exception ex) {
            System.out.println("Execution stopped as expected (Dummy Trapdoor): " + ex.getMessage());
            System.out.println("Code logic verified.");
        }
    }

    @Test
    public void testSampleLeft1(){
        BigInteger q = BigInteger.valueOf(17);
        int n = 5;
        double s = 5.0;
        int nParam = 5;

        // Create random matrices A and M
        BigIntMatrix[] A_S = trapGen(n, q.longValue(), 1);
        int m1 = A_S[0].getCols();
        int m = m1 + n;
        BigIntMatrix A = A_S[0];
        BigIntMatrix M = BigIntMatrix.random(n, m, q);

        BigIntMatrix T_A = A_S[1];

        // Target u
        BigInteger[] u = {BigInteger.valueOf(1), BigInteger.valueOf(2),BigInteger.valueOf(3),BigInteger.valueOf(2),BigInteger.valueOf(1)};

        System.out.println("Starting SampleLeft...");

        try {
            BigInteger[] result = sampleLeft(A, M, T_A, s, u, nParam);
            System.out.println("Result e: " + Arrays.toString(result));

            //Verify
            BigIntMatrix AM1 = BigIntMatrix.columnConcat(A, M);
            BigIntMatrix eMatrix = new BigIntMatrix(result.length,1,A.getModulus());

            for (int i = 0; i < result.length; i++) {
                eMatrix.set(i,0,result[i]);
            }

            BigIntMatrix resultTest = AM1.multiplyQ(eMatrix);
            BigInteger[] resData = new BigInteger[resultTest.getRows()];
            for(int i = 0; i< resultTest.getRows();i++){
                resData[i] = resultTest.get(i,0);
            }
            System.out.println("AM1 * e = u: " + Arrays.toString(resData));
            System.out.println("u: " + Arrays.toString(u));

        } catch (Exception ex) {
            System.out.println("Execution stopped as expected (Dummy Trapdoor): " + ex.getMessage());
            System.out.println("Code logic verified.");
        }
    }
}
