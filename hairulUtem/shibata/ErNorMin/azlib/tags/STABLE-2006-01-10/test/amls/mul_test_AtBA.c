#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include "mathematics.h"
#include "fem.h"
#include "base.h"



int main(void);

int main(void)
{
	int ii1, ii2, ii3, ii4, ii5, ii6;
	double val; 
	double **A;
	double **B;
	double **AtBA;
	double **AtB;
	double (**AtM)[3][3];
	double (**At)[3][3];
	double (**M)[3][3];
	double (**AtMA)[3][3];

	RC_TRY_MAIN( allocate2D(12, 12, &A) );
	RC_TRY_MAIN( allocate2D(12, 12, &B) );
	RC_TRY_MAIN( allocate2D(12, 12, &AtBA) );
	RC_TRY_MAIN( allocate2D(12, 12, &AtB) );
	RC_TRY_MAIN( allocate2D33(4, 4, &At) );
	RC_TRY_MAIN( allocate2D33(4, 4, &M) );
	RC_TRY_MAIN( allocate2D33(4, 4, &AtMA) );
	RC_TRY_MAIN( allocate2D33(4, 4, &AtM) );

	val = 0.0;
	for(ii1=0; ii1<12; ii1++){
		for(ii2=0; ii2<12; ii2++){
			val += 12.0;
			A[ii1][ii2] = val;
		}
	}

	val = 100.0;
	for(ii1=0; ii1<12; ii1++){
		for(ii2=0; ii2<ii1; ii2++){
			val += 20.0;
			B[ii1][ii2] = B[ii2][ii1] = val;
		}
		B[ii1][ii1] = val * 1.2;
	}

//	mul_matrix_AtB(16, 12, 12, A, C, AtC);
//	print_matrix(stderr, 6, 6, AtC);
	
	/* AtBA */
	for(ii1=0; ii1<12; ii1++){
		for(ii2=0; ii2<12; ii2++){
			for(ii3=0; ii3<12; ii3++){
				for(ii4=0; ii4<12; ii4++){
					AtBA[ii1][ii2] += A[ii3][ii1] * B[ii3][ii4] * A[ii4][ii2];
				}
			}
		}
	}
	/* AtB */
	for(ii1=0; ii1<12; ii1++){
		for(ii2=0; ii2<12; ii2++){
			for(ii3=0; ii3<12; ii3++){
				AtB[ii1][ii2] += A[ii3][ii1] * B[ii3][ii2];
			}
		}
	}

	/* B -> M */
	for(ii1=0; ii1<4; ii1++){
		for(ii2=0; ii2<=ii1; ii2++){
			for(ii3=0; ii3<3; ii3++){
				for(ii4=0; ii4<3; ii4++){
					M[ii1][ii2][ii3][ii4] = B[ii1*3+ii3][ii2*3+ii4];
				}
			}
		}
	}

	/* A -> At */
	for(ii1=0; ii1<4; ii1++){
		for(ii2=0; ii2<4; ii2++){
			for(ii3=0; ii3<3; ii3++){
				for(ii4=0; ii4<3; ii4++){
					At[ii1][ii2][ii3][ii4] = A[ii2*3+ii4][ii1*3+ii3];
				}
			}
		}
	}


#if 1
	for(ii1=0; ii1<4; ii1++){
		for(ii2=0; ii2<4; ii2++){
			for(ii3=0; ii3<ii2; ii3++){
				for(ii4=0; ii4<3; ii4++){
					for(ii5=0; ii5<3; ii5++){
						for(ii6=0; ii6<3; ii6++){
							AtM[ii1][ii2][ii5][ii4] += M[ii2][ii3][ii4][ii6]
							                        * At[ii1][ii3][ii5][ii6];
							AtM[ii1][ii3][ii5][ii4] += M[ii2][ii3][ii6][ii4]
							                        * At[ii1][ii2][ii5][ii6];
						}							
					}
				}
			}
//			for(ii3=ii2; ii3<ii2+1; ii3++){
				for(ii4=0; ii4<3; ii4++){
					for(ii5=0; ii5<3; ii5++){
						for(ii6=0; ii6<3; ii6++){
							AtM[ii1][ii2][ii5][ii4] += M[ii2][ii2][ii4][ii6]
							                        * At[ii1][ii2][ii5][ii6];
						}
					}
				}
//			}
		}
	}
	for(ii1=0; ii1<4; ii1++){
		for(ii2=0; ii2<4; ii2++){
			for(ii3=0; ii3<4; ii3++){
				for(ii4=0; ii4<3; ii4++){
					for(ii5=0; ii5<3; ii5++){
						for(ii6=0; ii6<3; ii6++){
							AtMA[ii1][ii2][ii4][ii5] += AtM[ii1][ii3][ii4][ii6]
							                          * At[ii2][ii3][ii5][ii6];
						}
					}
				}
			}
		}
	}

#endif
	/* Z_row */
	for(ii1=0; ii1<4; ii1++){
#if 0
		for(ii2=0; ii2<4; ii2++){
			for(ii3=0; ii3<ii2; ii3++){
				for(ii4=0; ii4<3; ii4++){
					for(ii5=0; ii5<3; ii5++){
						for(ii6=0; ii6<3; ii6++){
							AtM[ii2][ii5][ii4] += M[ii2][ii3][ii4][ii6]
							                    * At[ii1][ii3][ii5][ii6];
							AtM[ii3][ii4][ii5] += M[ii2][ii3][ii6][ii4]
							                    * At[ii1][ii2][ii5][ii6];
						}							
					}
				}
			}
			for(ii3=0; ii3<3; ii3++){
				for(ii4=0; ii4<3; ii4++){
					for(ii5=0; ii5<3; ii5++){
						AtM[ii2][ii3][ii4] += M[ii2][ii2][ii3][ii5]
						                    * At[ii1][ii2][ii4][ii5];
					}						
				}
			}
		}
#endif
#if 0
		for(ii2=0; ii2<4; ii2++){
			for(ii3=0; ii3<4; ii3++){
				for(ii4=0; ii4<3; ii4++){
					for(ii5=0; ii5<3; ii5++){
						for(ii6=0; ii6<3; ii6++){
							AtMA[ii1][ii2][ii4][ii5] +=
							AtM[ii3][ii4][ii6] * At[ii2][ii3][ii5][ii6];
						}							
					}
				}
			}
		}
#endif
	}



#if 0
	for(ii1=0; ii1<4; ii1++){
		for(ii2=0; ii2<4; ii2++){
			for(ii3=0; ii3<3; ii3++){
				for(ii4=0; ii4<3; ii4++){
					if(AtM[ii1][ii2][ii3][ii4] != AtB[ii1*3+ii3][ii2*3+ii4]){
					fprintf(stderr, "AtM[%d][%d][%d][%d] = %e\n",
					ii1,ii2,ii3,ii4, AtM[ii1][ii2][ii3][ii4]);
					fprintf(stderr,"AtB[%d][%d] = %e\n",ii1*3+ii3,ii2*3+ii4,
					AtB[ii1*3+ii3][ii2*3+ii4]);
					}
				}
			}
		}
	}
#endif
#if 1
	for(ii1=0; ii1<4; ii1++){
		for(ii2=0; ii2<4; ii2++){
			for(ii3=0; ii3<3; ii3++){
				for(ii4=0; ii4<3; ii4++){
					if(AtMA[ii1][ii2][ii3][ii4] != AtBA[ii1*3+ii3][ii2*3+ii4]){
					fprintf(stderr, "AtMA[%d][%d][%d][%d] = %e\n",
					ii1,ii2,ii3,ii4, AtMA[ii1][ii2][ii3][ii4]);
					fprintf(stderr,"AtBA[%d][%d] = %e\n",ii1*3+ii3,ii2*3+ii4,
					AtBA[ii1*3+ii3][ii2*3+ii4]);
					}
				}
			}
		}
	}
#endif





	RC_TRY_MAIN( free2D(12, 12, &A) );
	RC_TRY_MAIN( free2D(12, 12, &B) );
	RC_TRY_MAIN( free2D(12, 12, &AtBA) );
	RC_TRY_MAIN( free2D33(4, 4, &At) );
	RC_TRY_MAIN( free2D33(4, 4, &M) );
	RC_TRY_MAIN( free2D33(4, 4, &AtMA) );
	RC_TRY_MAIN( free2D33(4, 4, &AtM) );

	return(EXIT_SUCCESS);
}
