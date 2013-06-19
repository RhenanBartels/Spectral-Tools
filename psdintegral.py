 # _*_ coding: utf-8 _*_

 def psdint(F, Pxx, vlf_b=(0.003, 0.04),  lf_b=(0.04, 0.15), hf_b=(0.15, 0.4)):
     df = F[1] - F[0]
     vlf = trapz(Pxx[(F[F>=vlf_b[0]]) & (F[F <= vlf_b[1]])) * df
     lf = trapz(Pxx[(F[F>=lf_b[0]]) & (F[F <= lf_b[1]])) * df
     hf = trapz(Pxx[(F[F>=hf_b[0]]) & (F[F <= hf_b[1]])) * df
     TP = trapz(Pxx[(F[F>=vlf_b[0]]) & (F[F <= hf_b[1]])) * df
     lfnu = (lf / (TP - vlf)) * 100
     hfnu = (hf / (TP - vlf)) * 100
     return [vlf, lf, hf, TP, lfnu, hfnu]

