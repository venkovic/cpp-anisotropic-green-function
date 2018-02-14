int Binom(int n, int k) {
	if (k==0 || k==n) {
		return 1;
	} 
  return  Binom(n-1, k-1) + Binom(n-1, k);
}
