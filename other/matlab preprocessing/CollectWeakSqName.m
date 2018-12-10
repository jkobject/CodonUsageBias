SpeciesName={'sah6','saeu','saku','sp','so','sj'};

Name1=SequenceName1(Sid{1,1});
Name2=SequenceName2(Sid{1,2});
Name3=SequenceName3(Sid{1,3});
Name4=SequenceName4(Sid{1,4});
Name5=SequenceName5(Sid{1,5});
Name6=SequenceName6(Sid{1,6});

% setHeading(tbl, 'Asp Slection Pressure Analysis');

WeakSqAsp1 = table(Name1);
writetable(WeakSqAsp1);
WeakSqAsp2 = table(Name2);
writetable(WeakSqAsp2);
WeakSqAsp3 = table(Name3);
writetable(WeakSqAsp3);
WeakSqAsp4 = table(Name4);
writetable(WeakSqAsp4);
WeakSqAsp5 = table(Name5);
writetable(WeakSqAsp5);
WeakSqAsp6 = table(Name6);
writetable(WeakSqAsp6);

