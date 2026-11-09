function [LB0, LB] = safe_bound_error(blk, At, C, b, y, Z2, Bt, y2, l, max_eig)
  
  % Valid lower bounds through error bounds

  Aty = sdpnalAtyfun(blk, At, y);
  Znew = ops(C, '-', Aty);
  if ~isempty(Bt)
    Bty = sdpnalAtyfun(blk, Bt, y2);
    Znew = ops(Znew, '-', Bty);
  end
  if ~isempty(Z2)
     Znew = ops(Znew, '-', Z2); 
  end

  eigtmp = eig(full(Znew{1}));
  idx = find(eigtmp < -1e-8); 
  Xbar = max_eig;
  if isempty(l) || isempty(y2)
      LB0 = b'*y;
  else
      LB0 = b'*y + l'*y2;
  end
  pert = Xbar * sum(eigtmp(idx));
  LB = LB0 + pert; 
  
end