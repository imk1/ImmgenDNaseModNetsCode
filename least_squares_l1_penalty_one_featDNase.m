function [w0, w] = least_squares_l1_penalty_one_featDNase ( feat, featT, currentModMean, l1, l2, w )
% Use least squares with l1 and l2 penalties to learn the weights for the
% regulators

[nfeatures, nsamples] = size(feat);

w0 = mean(currentModMean);
currentModMean = currentModMean - w0;


a = sparse(1,nfeatures);
if(~exist('w'))
   % Initialize the weights to 0
   w = zeros(1,nfeatures);
else
   a( find( w~=0 ) ) = 1;
end
if(~exist('l2'))
   % Initialize the l2 penalty to 0
   l2 = zeros(nfeatures,1);
end
if( length(l1)==1 )
   % Doing flat Lirnet, so make a vector of the same repeated l1 penalties
   l1 = l1 * ones(nfeatures,1);
end
if( length(l2)==1 )
   % One universal l2 penalty, so make a vector of the same repeated l2
   % penalties
   l2 = l2 * ones(nfeatures,1);
end


alphaa = 0.3;
betaa = 0.8;

maxiters = 300;

for iter = 1 : maxiters
   % Iterate through features, considering the feature with the highest
   % current correlation on each iteration
    
   z = w * feat - currentModMean;

   preobj = sum( sum( z.^2 ) ) + abs(w)*l1 + (w.*w) * l2; % Compute the current objective

   grad = 2 * z * featT + 2 * l2' .* w; 

   abs_grad = abs(grad);
   abs_grad( find( a==1 ) ) = 0;

   [mx, mxindx] = max(abs_grad-l1');

   if( mx>0 )
      a(mxindx) = 1;w(mxindx) = 1.0e-5*sign(grad(mxindx));
   end

   active = find( a==1 );
   
   grad(active) = grad(active) + l1(active)' .* sign(w(active));

   grad2 = sum( sum( grad(active).^2 ) );

   t = 1;

   for biter = 1 : 300
   
      neww = w(active) - t * grad(active);

      newz = neww * feat(active,:) - currentModMean;

      postobj = sum( sum( newz.^2 ) ) + abs(neww)*l1(active) + (neww.*neww)*l2(active);

      if( postobj < preobj - alphaa * t * grad2 )
         break;
      else
         t = betaa * t;
      end

   end

   if( preobj > postobj )
      w(active) = neww;
   end

   if( abs(preobj-postobj) < 0.001 )
      break;
   end
end
