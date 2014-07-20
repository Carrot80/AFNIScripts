% compare Time intervalls:

Broca=[]
TimeInt600ms=[1 2 5 7 9 11 13 14 15 19 20 21];
TimeInt470ms=[3 6 8 10 12 16 17 18 22 23 24];

[p,h,stats]=signrank(mean(Broca(TimeInt600ms)), mean(Broca(TimeInt470ms)))
[p,h,stats]=signrank(mean(Wernicke(TimeInt600ms)), mean(Wernicke(TimeInt470ms)))