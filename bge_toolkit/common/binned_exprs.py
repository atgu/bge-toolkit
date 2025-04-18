import hail as hl


def mac_bin(s: hl.StructExpression):
    """MAC bin function.

    Assumes the field "AC" exists in the input StructExpression.

    Categories: ['1', '2-5', '6-10', '10+']
    """
    mac_bin = (hl.case()
               .when(s.AC <= 1, '1')
               .when(s.AC <= 5, '2-5')
               .when(s.AC <= 10, '6-10')
               .default('10+'))
    return mac_bin


def maf_bin(s: hl.StructExpression):
    """MAF bin function.

    Assumes the field "AF" exists in the input StructExpression.

    Categories: ['<1%', '1-2%', '2-5%', '>5%']
    """
    maf_bin = (hl.case()
               .when(s.AF <= 0.01, '<1%')
               .when(s.AF <= 0.02, '1-2%')
               .when(s.AF <= 0.05, '2-5%')
               .default('>5%'))
    return maf_bin


def gq_bin(s: hl.StructExpression):
    """GQ bin function.

    Assumes the field "GQ" exists in the input StructExpression.

    Categories: ['10', '20', '30', '40', '50', '60', '70', '80', '90']
    """
    x = (hl.case()
         .when(s.GQ <= 10, 10)
         .when(s.GQ <= 20, 20)
         .when(s.GQ <= 30, 30)
         .when(s.GQ <= 40, 40)
         .when(s.GQ <= 50, 50)
         .when(s.GQ <= 60, 60)
         .when(s.GQ <= 70, 70)
         .when(s.GQ <= 80, 80)
         .default(90))
    return x


def dp_bin(s: hl.StructExpression):
    """DP bin function.

    Assumes the field "DP" exists in the input StructExpression.

    Categories: ['10', '20', '30', '40', '50', '60', '70', '80', '90']
    """
    x = (hl.case()
         .when(s.DP <= 10, 10)
         .when(s.DP <= 20, 20)
         .when(s.DP <= 30, 30)
         .when(s.DP <= 40, 40)
         .when(s.DP <= 50, 50)
         .when(s.DP <= 60, 60)
         .when(s.DP <= 70, 70)
         .when(s.DP <= 80, 80)
         .when(s.DP <= 90, 90)
         .default(100))
    return x


def max_gp_bin(s: hl.StructExpression):
    """Max GP bin function.

    Assumes the field "MAX_GP" exists in the input StructExpression.

    Categories: ['0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '0.95', '1.0']
    """
    x = (hl.case()
         .when(s.MAX_GP <= 0.1, 0.1)
         .when(s.MAX_GP <= 0.2, 0.2)
         .when(s.MAX_GP <= 0.3, 0.3)
         .when(s.MAX_GP <= 0.4, 0.4)
         .when(s.MAX_GP <= 0.5, 0.5)
         .when(s.MAX_GP <= 0.6, 0.6)
         .when(s.MAX_GP <= 0.7, 0.7)
         .when(s.MAX_GP <= 0.8, 0.8)
         .when(s.MAX_GP <= 0.9, 0.9)
         .when(s.MAX_GP <= 0.95, 0.95)
         .default(1.0))
    return x


def qual_approx_bin(s: hl.StructExpression):
    """Approx Qual bin function.

    Assumes the field "QUALapprox" exists in the input StructExpression.

    Categories: ['10', '20', '30', '40', '50', '60', '70', '80', '90']
    """
    x = (hl.case()
         .when(s.QUALapprox <= 10, 10)
         .when(s.QUALapprox <= 20, 20)
         .when(s.QUALapprox <= 30, 30)
         .when(s.QUALapprox <= 40, 40)
         .when(s.QUALapprox <= 50, 50)
         .when(s.QUALapprox <= 60, 60)
         .when(s.QUALapprox <= 70, 70)
         .when(s.QUALapprox <= 80, 80)
         .default(90))
    return x


def info_score_bin(s: hl.StructExpression):
    """INFO score bin function.

    Assumes the field "INFO" exists in the input StructExpression.

    Categories: ['0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '0.95', '1.0']
    """
    x = (hl.case()
         .when(s.INFO <= 0.1, 0.1)
         .when(s.INFO <= 0.2, 0.2)
         .when(s.INFO <= 0.3, 0.3)
         .when(s.INFO <= 0.4, 0.4)
         .when(s.INFO <= 0.5, 0.5)
         .when(s.INFO <= 0.6, 0.6)
         .when(s.INFO <= 0.7, 0.7)
         .when(s.INFO <= 0.8, 0.8)
         .when(s.INFO <= 0.9, 0.9)
         .when(s.INFO <= 0.95, 0.95)
         .default(1.0))
    return x
