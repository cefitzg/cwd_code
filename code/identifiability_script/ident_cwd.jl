using StructuralIdentifiability


cwd_model= @ODEmodel(

    S'(t) = a -S*(b*J+g*E+m),

    J'(t) = S*(b*J+g*E)-J*(m+mu),

    E'(t) = ep*J-tau*E,

    C'(t) = mu*J,

    y1(t) = C(t)


    )

print(assess_identifiability(cwd_model))