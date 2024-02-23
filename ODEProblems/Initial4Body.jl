function Initial4Body(T=Float64)


   Gm0 = parse(BigFloat,"1")
   Gm1 = parse(BigFloat,"1e-3")
   Gm2 = parse(BigFloat,"1e-2")
   Gm3 = parse(BigFloat,"1e-2")

   Gm = [Gm0, Gm1, Gm2, Gm3]

   q = [parse(BigFloat,"0.016194591905998788"), 
        parse(BigFloat,"-0.0004991194886550273"),
        parse(BigFloat,"0.0"),

        parse(BigFloat,"-0.5790148460298099"),
        parse(BigFloat,"0.8444732531667211"),
        parse(BigFloat,"0.0"),

        parse(BigFloat,"-0.5503975770982608"),
        parse(BigFloat,"-0.8513300618402505"),
        parse(BigFloat,"0.0"),

        parse(BigFloat,"-1.011160128898637"),
        parse(BigFloat,"0.8167946853890811"),
        parse(BigFloat,"0.0")
      ]

   v = [parse(BigFloat,"0.014430763618203299"), 
      parse(BigFloat,"0.0018395077015696387"),
      parse(BigFloat,"0.0"),

      parse(BigFloat,"-0.7548037657105789"),
      parse(BigFloat,"-0.5528187854004273"),
      parse(BigFloat,"0.0"),

      parse(BigFloat,"-0.8169192926401513"),
      parse(BigFloat,"0.5570022821589271"),
      parse(BigFloat,"0.0"),

      parse(BigFloat,"-0.5506766926091207"),
      parse(BigFloat,"-0.6856711737758482"),
      parse(BigFloat,"0.0")
    ]


      q0 = reshape(q,3,:)
      v0 = reshape(v,3,:)

#      q0bar = [sum(Gm .* q0[j,:])/sum(Gm) for j in 1:3]
#      v0bar = [sum(Gm .* v0[j,:])/sum(Gm) for j in 1:3]

#      q0 = q0 .- q0bar
#      v0 = v0 .- v0bar

      N = length(Gm)    

      u0 = Array{T}(undef,3,N,2)
      u0[:,:,2] .= v0[:,1:N]
      u0[:,:,1] .= q0[:,1:N]

      bodylist = ["body1" "body2" "body3" "body4"]

      return u0, Gm, bodylist

end

