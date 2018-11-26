/** 
	hopfield.pl 
	
	Τρέχει και σε SWI, και σε ECLiPse.
*/

/** 
	hopfield calculation 
	
	Αφού βρει τα M, N, υπολογίζει τον μοναδιάιο πίνακα, και το M * In. Στη σύνέχεια, υπολογίζει (το υπόλοιπο κομμάτι),
	τους πολλαπλασιασμούς διανυσμάτων και τα αθροίσματα αυτών (matlist, mySum) και επιστρέφει τη διαφορά αυτών των 
	δύο κομματιών στη μεταβλητή Α. Τέλος χτίζεται η κατάλληλη λίστα, και το αποτέλεσμα επιστρέφεται στη W.
*/
hopfield(Ls, W) :- length(Ls, M), first(Ls, F), length(F, N), id(N, I), mulLs(M, I, Id), matlist(Ls, L), mySum(L, NewList), flatten(Id, NewId), diff2(NewList, NewId, A), build(N, A, W).

/** Επιστρέφει σε μία (flat) λίστα το αποτέλεσμα όλων των πολλαπλασιασμών των διανυσμάτων που δίνονται στην (λίστα λιστών) Ls. */ 
matlist(Ls, L) :- matlist(Ls, [], L).
matlist([], Acc, Acc).
matlist([H|T], Acc, L) :- matmul(H, M), flatten(M, K), append(Acc, [K], NewList), matlist(T, NewList, L).


/** 
	Identity matrix

	Υπολογίζει τον μοναδιαίο πίνακα για το συγκεκριμένο μέγεθος N που δίνουμε. Επιστρέφει το αποτέλεσμα
	σε μία λίστα λιστών. Το κατηγόρημα zeros γεμίζει κάθε υπολίστα με τον δοθέν αριθμό μηδενικών και στο id, 
	αντικαθίσταται με 1 η κατάλληλη θέση κάθε φορά.
*/
id(N, L) :- id(N, 0, [], L).
id(N, N, Acc, Acc).
id(N, Temp, Acc, L) :- Temp =\= N, Count is Temp +1, zeros(N, K), rep(K, Temp, 1, NewList), append(Acc, [NewList], List), id(N, Count, List, L).

zeros(N, L) :- zeros(N, 0, [], L).
zeros(N, N, Acc, Acc).
zeros(N, Temp, Acc, L) :- Temp =\= N, Count is Temp +1, append(Acc, [0], K), zeros(N, Count, K, L).


/** 
	Matrix multiplication

	Πολλαπλασιάζει το διάνυσμα (λίστα) του δίνουμε, με το ανάστροφό του, και επιστρέφει το αποτέλεσμα
	σε μία λίστα (λιστών). Με τη βοήθεια των findelem, muln, πολλαπλασιάζει κάθε στοιχείο με κάθε στοιχείο
	της δοθείσας λίστας.
*/
matmul(Ls, L) :- matmul(Ls, 0, [], L).
matmul(Ls, Temp, Acc, Acc) :- length(Ls, J), Temp > J-1.
matmul(Ls, Temp, Acc, L) :- findelem(Ls, Temp, F), Count is Temp + 1, muln(F, Ls, K), append(Acc, [K], Z), matmul(Ls, Count, Z, L).  


/** 

	Number with list multiplication 

	muln : Πολλαπλασιάζει έναν αριθμό N, με κάθε στοιχείο της λίστας Ls.
	mulLs : Πολλαπλασιάζει έναν αριθμό N, με κάθε στοιχείο κάθε λίστας της λίστας (λιστών) Ls.

*/
muln(N, Ls, L) :- muln(N, Ls, [], L).
muln(_, [], Acc, Acc).
muln(X, [H|T], Acc, K) :- Temp is X * H, append(Acc, [Temp], B), muln(X, T, B, K).
 
mulLs(N, Ls, L) :- mulLs(N, Ls, [], L).
mulLs(_, [], Acc, Acc).
mulLs(N, Ls, Acc, L) :- first(Ls, H), muln(N, H, F), append(Acc, [F], P), ini(Ls, I), mulLs(N, I, P, L).


/**
	 General scope predicates 
*/

/** Επιστρέφει το στοιχείο της δοθεισας λιστας στη θέση I. */
findelem([H|_], 0, H).
findelem([_|T], I, X):- Index is I - 1, findelem(T, Index, X).

/** Αντικαθιστά τη θέση I, με το στοιχείο X, στη λίστα με κεφαλή το H, και ουρά το T. */
rep([_|T], 0, NewHead, [NewHead|T]). 
rep([H|T], I, X, [H|NewTail]):- Index is I - 1, rep(T, Index, X, NewTail).

/** Χρησιμοποιούνται για την επιστροφή του head (η first), και tail (η ini), σε ορισμένες περιπτώσεις που χρειάζεται. */
first([H|_], H).
ini([_|T], T).

/** Διαφορά στοιχείων (στοιχείο με στοιχείο, όχι συνολική) δύο λιστών. Επιστρέφει το αποτέλεσμα σε μία νέα λίστα. */
diff2(L1, L2, L) :- diff2(L1, L2, [], L).
diff2([], [], Acc, Acc).
diff2([H1|T1], [H2|T2], Acc, L) :- Head is H1-H2, append(Acc, [Head], NewList), diff2(T1, T2, NewList, L).

/** Άθροισμα στοιχείων (στοιχείο με στοιχείο, όχι συνολικό) δύο λιστών. Επιστρέφει το αποτέλεσμα σε μία νέα λίστα. */
sum2(L1, L2, L) :- sum2(L1, L2, [], L).
sum2([], [], Acc, Acc).
sum2([H1|T1], [H2|T2], Acc, L) :- Head is H1+H2, append(Acc, [Head], NewList), sum2(T1, T2, NewList, L).

/** 
	"Σπάει" τη δοθείσα λίστα σε N υπολίστες, και χωρίζει σε N στοιχεία τη καθεμία. Επιστρέφει σε μία λίστα το 
	 αποτέλεσμα, οπότε στο τέλος έχουμε μία λίστα λιστών. Αν δεν βρει τον κατάλληλο αριθμό στοιχείων για 
	 ισομοιρασμό στις υπολίστες, επιστρέφει No.

	 π.χ. ?-build(3, [4,5,6,3,45,3,1,2,1], M).
			M = [[4, 5, 6], [3, 45, 3], [1, 2, 1]]; No.

		  ?-build(2, [4,5,6,3,45,3,1,2,1], M).
		    No.
*/
build(N, Ls, W) :- build(N, 0, Ls, [], [], W).
build(_, 0, [], [], Acc2, Acc2). 
build(N, Temp, [H|T], Acc1, Acc2, W) :- Temp =\= N, Count is Temp +1, append(Acc1, [H], K), buildme(N, Count, T, K, Acc2, W).

buildme(N, Count, T, K, Acc2, W) :- Count =\= N, build(N, Count, T, K, Acc2, W).
buildme(N, N, T, K, Acc2, W) :- append(Acc2, [K], L), build(N, 0, T, [], L, W).

/** Άθροισμα στοιχείων (στοιχείο με στοιχείο, όχι συνολικό) κάθε υπολίστας της λίστας (λιστών) που δίνουμε. Επιστρέφει το αποτέλεσμα σε μία νέα λίστα. */
mySum([H|[]], H).
mySum([H|T], L) :- first(T, F), sum2(H, F, TempList), ini(T, I), mySum([TempList|I], L).  