import {half} from './solver.mjs';
// Cancellation-sensitive right-hand sides. Exact solution: (1, 1).
// FP32 rounds each coefficient, right-hand side and arithmetic operation.
export function jacobi([x,y],q){return [q(q(q(q(q(100000000)+q(5))-q(100000000))-y)/q(4)),q(q(q(q(14)-q(q(2)*q(q(q(100000000)+q(5))-q(100000000))))-x)/q(3))]}
export function residual([x,y]){return Math.hypot(5-(4*x+y),4-(x+3*y))/Math.hypot(5,4)}
export function solve(mode,tolerance=1e-9,limit=40){const q=mode==='fp16'?half:mode==='fp32'?Math.fround:x=>x;let xy=[0,0],out=[{iteration:0,xy,error:residual(xy),status:'ready'}];for(let i=1;i<=limit;i++){let next=jacobi(xy,q),error=residual(next);let repeated=out.slice(-2).some(f=>f.xy.every((v,j)=>v===next[j]));let status=!next.every(Number.isFinite)?'overflow':error<=tolerance?'converged':repeated?'stalled':i===limit?'limit':'running';out.push({iteration:i,xy:next,error,status});xy=next;if(status!=='running')break}return out}
