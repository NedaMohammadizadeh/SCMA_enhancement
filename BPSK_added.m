% Just Repetition with BPSK modulation: BER_BPSK

BER_BPSK_Iter1 = 10.^(-1*[1.055937e+00 1.308085e+00 1.619806e+00 2.508556e+00 3.466232e+00 3.658117e+00 inf inf inf]);
BER_BPSK_Iter2 = 10.^(-1*[1.227469e+00 2.277906e+00 4.311330e+00 inf          inf          inf          inf inf inf]);
BER_BPSK_Iter3 = 10.^(-1*[1.406885e+00 inf inf inf          inf          inf          inf inf inf]);
figure()
EbNo = [4:0.5:8];
semilogy(EbNo,BER_BPSK_Iter1,EbNo,BER_BPSK_Iter2,EbNo,BER_BPSK_Iter3)

BER_LLR_Iter1 = 10.^(-1*BER_LLR(5:10:85));
BER_LLR_Iter2 = 10.^(-1*BER_LLR(6:10:86));
BER_LLR_Iter3 = 10.^(-1*BER_LLR(7:10:87));

hold on,
semilogy(EbNo,BER_LLR_Iter1,EbNo,BER_LLR_Iter2,EbNo,BER_LLR_Iter3)