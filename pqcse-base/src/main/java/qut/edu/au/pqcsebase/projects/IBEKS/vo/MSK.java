package qut.edu.au.pqcsebase.projects.IBEKS.vo;

import lombok.*;
import qut.edu.au.pqcsebase.tools.BigIntMatrix;

/**
 * Master Secret Key class for IBEKS scheme.
 */
@Data@Setter@Getter@AllArgsConstructor@ToString
public class MSK {

    // m*m dimension matrix basis
    public BigIntMatrix msk;
}
