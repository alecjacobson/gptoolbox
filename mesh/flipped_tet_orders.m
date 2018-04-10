function flipped_order = flipped_tet_orders()
  % FLIPPED_TET_ORDERS
  %
  % flipped_order = flipped_tet_orders()
  %
  % Outputs:
  %   flipped_order  20 by 4 list of tet index orders that are flipped inside
  %     out (negative volume)
  %
  flipped_order = [ ...
      4 3 1 2
      4 2 3 1
      4 1 2 3
      3 4 2 1
      3 2 1 4
      3 1 4 2
      2 4 1 3
      2 3 4 1
      2 1 3 4
      1 4 3 2
      1 3 2 4
      1 2 4 3];
end
