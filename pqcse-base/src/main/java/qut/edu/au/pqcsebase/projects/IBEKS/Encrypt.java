package qut.edu.au.pqcsebase.projects.IBEKS;

import qut.edu.au.pqcsebase.projects.IBEKS.vo.CT;
import qut.edu.au.pqcsebase.tools.BigIntMatrix;

import java.math.BigInteger;
import java.security.SecureRandom;
import java.util.Arrays;

import static qut.edu.au.pqcsebase.algorithms.ErrorDistribution.generateError;
import static qut.edu.au.pqcsebase.projects.IBEKS.Params.N;
import static qut.edu.au.pqcsebase.projects.IBEKS.Params.n;
import static qut.edu.au.pqcsebase.tools.BigIntMatrix.*;
import static qut.edu.au.pqcsebase.tools.VectorUtils.*;

public class Encrypt {

    public static CT encrypt(Setup setup, String identity, String keyword) {

        BigIntMatrix A = setup.pp.A;

        BigIntMatrix R_id = setup.hash1(identity);
        BigIntMatrix A_id = A.multiplyQ(inverseQ(R_id));

        //keyword
        BigInteger x_w = setup.hash2(keyword);

        BigInteger[] y0 = new BigInteger[(int)n];
        for (int i = 0; i < (int)N + 1; i++) {
            y0[i] = x_w.modPow(BigInteger.valueOf(i), BigInteger.valueOf(Params.q));
        }
        System.out.println("y0: " + Arrays.toString(y0));

        BigInteger[] y1 = setup.pp.y1;

        //[y0|y1]
        for(int i = 0; i< y1.length; i++) {
            y0[(int)N + 1 + i] = y1[i];
        }

        //generate v and s
        BigInteger[] v = randomVector((int)n, BigInteger.valueOf(Params.q));
        BigInteger[] s = randomVector((int)n, BigInteger.valueOf(Params.q));

        //generate R
        int m = A.getColumnDimension();
        BigIntMatrix R = new BigIntMatrix(m,m, BigInteger.valueOf(Params.q));
        SecureRandom secureRandom = new SecureRandom();
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                R.set(i,j, secureRandom.nextBoolean()? BigInteger.ONE : BigInteger.valueOf(-1));
            }
        }
        System.out.println("R size: " + R.getRowDimension() + "x" + R.getColumnDimension());
        System.out.println("R:\n" + R);

        //generate errors
        BigInteger z_u = generateError(BigInteger.valueOf(10007));
        BigInteger[] z_id = new BigInteger[m];
        for (int i = 0; i < m; i++) {
            z_id[i] = generateError(BigInteger.valueOf(10007));
        }

        BigInteger[] z_w = multiplyMatrixVector(R,z_id);
        System.out.println("z_id.length: " + z_id.length);
        System.out.println("z_id: " + Arrays.toString(z_id));
        System.out.println("z_w: " + Arrays.toString(z_w));

        //compute ciphertext
        BigInteger c_u = innerProduct(setup.pp.u, s).add(z_u).mod(BigInteger.valueOf(Params.q));
        BigInteger[] c_id = addVectorsMod(multiplyMatrixVector(A_id.transpose(), s), z_id, BigInteger.valueOf(Params.q));

        BigInteger vy = innerProduct(v, y0).mod(BigInteger.valueOf(Params.q));
        System.out.println("vy: " + vy);
        BigIntMatrix B_vyA = new BigIntMatrix((int)n,m,BigInteger.valueOf(Params.q));
        System.out.println("B_vyA size: " + B_vyA.getRowDimension() + "x" + B_vyA.getColumnDimension());

        for(int i = 0; i < (int)n; i++){
            for(int j = 0; j < m; j++){
                BigInteger val = A_id.get(i,j).multiply(vy);
                //System.out.println("val at (" + i + "," + j + "): " + val);
                B_vyA.set(i,j, setup.pp.B.get(i,j).add(val).mod(BigInteger.valueOf(Params.q)));
                //System.out.println("B_vyA at (" + i + "," + j + "): " + B_vyA.get(i,j));
            }
        }
        System.out.println("B_vyA:\n" + B_vyA);

        BigInteger[] c_w = addVectorsMod(multiplyMatrixVector(B_vyA.transpose(), s), z_w, BigInteger.valueOf(Params.q));

        return new CT(c_u,c_id,c_w);
    }

}
