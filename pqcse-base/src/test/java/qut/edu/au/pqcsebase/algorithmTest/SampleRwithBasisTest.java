package qut.edu.au.pqcsebase.algorithmTest;

import org.junit.Test;
import qut.edu.au.pqcsebase.tools.BigIntMatrix;

import static java.lang.Math.log;
import static java.lang.Math.sqrt;
import static qut.edu.au.pqcsebase.algorithms.SampleRwithBasis.sampleRwithBasis;
import static qut.edu.au.pqcsebase.algorithms.TrapGen.trapGen;
import static qut.edu.au.pqcsebase.tools.BigIntMatrix.inverseQ;

public class SampleRwithBasisTest {

    @Test
    public void testSampleRwithBasis() {
        long n = 5;
        long q = 17;
        int flag = 1;
        BigIntMatrix[] A_S = trapGen(n,q,flag);
        BigIntMatrix A = A_S[0];
        BigIntMatrix S = A_S[1];
        int m = S.getRowDimension();
        System.out.println("m :" + m);
        double sigmaR = sqrt(n * log(q)) * sqrt(log(m));
        double sigma = 10000;
        System.out.println("sigmaR: " + sigmaR);
        BigIntMatrix[] R_Tb = sampleRwithBasis(A, sigmaR);
        BigIntMatrix R = R_Tb[0];
        BigIntMatrix Tb = R_Tb[1];
        // B = B = AR^(-1) mod q
        BigIntMatrix B = A.multiplyQ(inverseQ(R));
        System.out.println("B: \n" + A.toString());
        System.out.println("R: \n" + R.toString());
        System.out.println("Tb: \n" + Tb.toString());
        BigIntMatrix result = B.multiplyQ(Tb);
        System.out.println("B * Tb mod q: \n" + result.toString());
        System.out.println(result.isZeroMatrix());
    }
}
