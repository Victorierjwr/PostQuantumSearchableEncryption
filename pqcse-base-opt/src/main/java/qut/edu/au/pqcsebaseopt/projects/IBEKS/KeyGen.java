package qut.edu.au.pqcsebaseopt.projects.IBEKS;


import qut.edu.au.pqcsebaseopt.projects.IBEKS.vo.KeyPair;
import qut.edu.au.pqcsebaseopt.projects.IBEKS.vo.MSK;
import qut.edu.au.pqcsebaseopt.tools.BigIntMatrix;


import static qut.edu.au.pqcsebaseopt.algorithms.BasisDel.basisDel;
import static qut.edu.au.pqcsebaseopt.projects.IBEKS.Params.n;
import static qut.edu.au.pqcsebaseopt.projects.IBEKS.Params.sigma;
import static qut.edu.au.pqcsebaseopt.tools.BigIntMatrix.inverseQ;

public class KeyGen {

    public static KeyPair keyGen(Setup setup, MSK msk, String identity) {

        BigIntMatrix A = setup.pp.A;
        BigIntMatrix TA = msk.msk;
        BigIntMatrix R_id = setup.hash1(identity);
        System.out.println("R_id: " + R_id);

        BigIntMatrix R_id_inv = inverseQ(R_id);
        System.out.println("R_id_inv: " + R_id_inv);

        BigIntMatrix A_id = A.multiplyQ(R_id_inv);

        BigIntMatrix T_id = basisDel(A,R_id,TA,sigma,n);

        return new KeyPair(A_id,T_id);
    }
}
