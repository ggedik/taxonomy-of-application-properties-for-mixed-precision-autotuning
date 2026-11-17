// Binary16, round to nearest with ties to even. Arithmetic is emulated in JS.
export function half(value) {
  if (!Number.isFinite(value) || value === 0) return value;
  const sign = Math.sign(value), a = Math.abs(value);
  if (a >= 65520) return sign * Infinity;
  const spacing = 2 ** Math.max(-24, Math.floor(Math.log2(a)) - 10);
  const scaled = a / spacing, lower = Math.floor(scaled), fraction = scaled - lower;
  const integer = fraction > 0.5 || (fraction === 0.5 && lower % 2) ? lower + 1 : lower;
  return sign * integer * spacing;
}
export const rounders = { fp64: x => x, fp32: Math.fround, fp16: half };
export function residual([x, y]) {
  return Math.hypot(1 - (4 * x + y), 2 - (x + 3 * y)) / Math.sqrt(5);
}
// Simultaneous Jacobi updates for 4x + y = 1, x + 3y = 2.
export function jacobi([x, y], q) {
  return [q(q(q(1) - y) / q(4)), q(q(q(2) - x) / q(3))];
}
export function solve(mode, tolerance = 1e-6, limit = 24) {
  const q = rounders[mode];
  if (!q) throw new Error('Unknown precision');
  let xy = [q(0), q(0)];
  const frames = [{ iteration: 0, xy, error: residual(xy), status: 'ready' }];
  for (let iteration = 1; iteration <= limit; iteration++) {
    const next = jacobi(xy, q), error = residual(next);
    const repeated = frames.slice(-2).some(f => f.xy.every((v, i) => v === next[i]));
    const status = error <= tolerance ? 'converged' : repeated ? 'stalled' : iteration === limit ? 'limit' : 'running';
    frames.push({ iteration, xy: next, error, status });
    xy = next;
    if (status !== 'running') break;
  }
  return frames;
}
