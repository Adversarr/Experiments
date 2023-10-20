clear;
N = 10;
NStep = 1000;
f = zeros(N, N, 9);
f(2, 2, :) = ones(1, 9);
figure;
hm = heatmap(compute_rho(f))
pause(2)
ci = int32([0 0; 1 0; 0 1; 1 1; -1 0; 0 -1; -1 -1; 1 -1; -1 1]);
for i = 1: NStep
  f = prop(f, N, ci);
  hm.ColorData = compute_rho(f);
  sum(f, "all")
  pause(2);
end

function y = compute_rho(f)
  y = sum(f, 3);
end

function fp = prop(f, N, ci)
  fp = zeros(N, N, 9);
  for i = 2: N-1
    for j = 2:N-1
      for k = 1:9
        fp(i+ci(k, 1), j+ci(k, 2), k) = fp(i+ci(k, 1), j+ci(k, 2), k) ...
          + f(i, j, k);
      end
    end
  end
  % deal with boundaries:
  for i = 2: N-1
    for k = [1, 2, 3, 4, 5, 9]
      fp(i + ci(k, 1), 1+ci(k, 2), k) = fp(i + ci(k, 1), 1+ci(k, 2), k) ...
        + f(i, 1, k);
    end
    for k = [1, 2, 5, 6, 7, 8]
      fp(i + ci(k, 1), N+ci(k, 2), k) = fp(i + ci(k, 1), N+ci(k, 2), k) ...
        + f(i, N, k);
    end
  end
  for j = 2:N-1
    for k = [1, 2, 3, 4, 6, 8]
      fp(1 + ci(k, 1), j+ci(k, 2), k) = fp(1 + ci(k, 1), j+ci(k, 2), k) ...
        + f(1, j, k);
    end
    for k = [1, 3, 5, 6, 7, 9]
      fp(N + ci(k, 1), j+ci(k, 2), k) = fp(N + ci(k, 1), j+ci(k, 2), k) ...
        + f(N, j, k);
    end
  end
  % i = j = 1
  for k = 1:4
    fp(1 + ci(k, 1), 1+ci(k, 2), k) = fp(1 + ci(k, 1), 1+ci(k, 2), k) ...
        + f(1, 1, k);
  end
  % i = j = N
  for k = [1 5 6 7]
    fp(N + ci(k, 1), N+ci(k, 2), k) = fp(N + ci(k, 1), N+ci(k, 2), k) ...
        + f(N, N, k);
  end
  % i = n, j = 1
  for k = [1 3 5 9]
    fp(N + ci(k, 1), 1+ci(k, 2), k) = fp(N + ci(k, 1), 1+ci(k, 2), k) ...
        + f(N, 1, k);
  end
  % j = n, i = 1
  for k = [1 2 6 8]
    fp(1 + ci(k, 1), N+ci(k, 2), k) = fp(1 + ci(k, 1), N+ci(k, 2), k) ...
        + f(1, N, k);
  end
end

