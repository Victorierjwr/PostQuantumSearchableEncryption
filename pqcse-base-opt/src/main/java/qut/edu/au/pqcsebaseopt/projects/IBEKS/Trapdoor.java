package qut.edu.au.pqcsebaseopt.projects.IBEKS;


import qut.edu.au.pqcsebaseopt.projects.IBEKS.vo.TW;
import qut.edu.au.pqcsebaseopt.tools.BigIntMatrix;


import java.math.BigInteger;
import java.security.SecureRandom;
import java.util.Arrays;

import static qut.edu.au.pqcsebaseopt.algorithms.SampleLeft.sampleLeft;
import static qut.edu.au.pqcsebaseopt.algorithms.SamplePre.samplePre;
import static qut.edu.au.pqcsebaseopt.projects.IBEKS.Params.*;
import static qut.edu.au.pqcsebaseopt.tools.BigIntMatrix.*;
import static qut.edu.au.pqcsebaseopt.tools.PolynomialUtils.getCoeffsFromRoots;


public class Trapdoor {

    public static TW trapdoor(Setup setup, String[] keywords, BigIntMatrix sk_id, String identity){
        BigIntMatrix A = setup.pp.A;

        BigIntMatrix R_id = setup.hash1(identity);
        BigIntMatrix A_id = A.multiplyQ(inverseQ(R_id));

        int l = keywords.length;
        BigInteger[] H_w = new BigInteger[(int)N];
        for(int i=0;i<l;i++){
            H_w[i] = setup.hash2(keywords[i]);
        }

        for (int i = l; i < N; i++) {
            SecureRandom secureRandom = new SecureRandom();
            H_w[i] = BigInteger.valueOf(secureRandom.nextInt()).mod(BigInteger.valueOf(Params.q));
        }
        System.out.println("H_w: " + Arrays.toString(H_w));

        BigInteger[] b0 = getCoeffsFromRoots(H_w, BigInteger.valueOf(Params.q));
        System.out.println("b0: " + Arrays.toString(b0));
        BigInteger[] b = new BigInteger[(int)n];
        System.arraycopy(b0, 0, b, 0, b0.length);
        for(int i=b0.length;i< n;i++){
            b[i] = setup.pp.b1[i-b0.length];
        }

        BigInteger[] e0 = samplePre(A_id, sk_id, sigma,b);
        System.out.println("e0: " + Arrays.toString(e0));
        System.out.println("A_id * e0 == b mod q: " + Arrays.toString(multiplyMatrixVector(A_id,e0)));
        System.out.println("b mod q: " + Arrays.toString(b));

        BigInteger[] Be0_tmp = multiplyMatrixVector(setup.pp.B,e0);
        BigIntMatrix Be0 = new BigIntMatrix(Be0_tmp.length,1,BigInteger.valueOf(Params.q));
        for(int i=0;i<Be0_tmp.length;i++){
            Be0.set(i,0,Be0_tmp[i]);
        }

        BigInteger[] e1 = sampleLeft(A_id,Be0,sk_id,sigma,setup.pp.u,(int)n);
        System.out.println("e1: " + Arrays.toString(e1));
        BigIntMatrix tmp = columnConcat(A_id,Be0);
        System.out.println("[A_id|Be0] * e1== u mod q: " + Arrays.toString(multiplyMatrixVector(tmp,e1)));
        System.out.println("u mod q: " + Arrays.toString(setup.pp.u));

        return new TW(e0,e1);

    }
}
