function data=linesearch_simple(funfcn, data, optim)
% Find a bracket of acceptable points
data = bracketingPhase_simple(funfcn, data, optim);
if (data.bracket_exitflag  == 2)
  % BracketingPhase found a bracket containing acceptable points; 
  % now find acceptable point within bracket
  data = sectioningPhase_simple(funfcn, data, optim);
  data.exitflag = data.section_exitflag; 
else
  % Already acceptable point found or MaxFunEvals reached
  data.exitflag = data.bracket_exitflag; 
end
