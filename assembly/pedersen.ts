import { BigInt } from './bigint';

// Convert a string representation of a number to a BigInt object
function BN(n: string): BigInt {
  return BigInt.fromString(n);
}

// https://docs.starkware.co/starkex/stark-curve.html

// Params: a, b
export const CURVE_A = BN('1');
export const CURVE_B = BN(
  '3141592653589793238462643383279502884197169399375105820974944592307816406665'
);

// Field over which we'll do calculations. Verify with:
// NOTE: there is no efficient sqrt for field (P%4==1)
export const CURVE_P = BN('2')
  .pow(251)
  .add(BN('17').mul(BN('2').pow(192)))
  .add(BN('1'));

// Field over which we'll do calculations. Verify with:
// NOTE: there is no efficient sqrt for field (P%4==1)
// P: BN('2').pow(251).add(BN('17').mul(BN('2').pow(192)).add(BN('1'))),
// Curve order, total count of valid points in the field. Verify with:
export const CURVE_N = BN(
  '3618502788666131213697322783095070105526743751716087489154079457884512865583'
);
export const CURVE_N_BITS = 252; // len(bin(N).replace('0b',''))

// Base point (x, y) aka generator point
export const CURVE_GX = BN(
  '874739451078007766457464989774322083649278607533249481151382481072868806602'
);
export const CURVE_GY = BN(
  '152666792071518830868575557812948353041420400780739481342941381225525861407'
);

/**
 * Jacobian Point works in 3d / jacobi coordinates: (x, y, z) ∋ (x=x/z², y=y/z³)
 * Default Point works in 2d / affine coordinates: (x, y)
 * We're doing calculations in jacobi, because its operations don't require costly inversion.
 */
class JacobianPoint {
  constructor(readonly x: BigInt, readonly y: BigInt, readonly z: BigInt) {}

  static readonly BASE: JacobianPoint = new JacobianPoint(
    CURVE_GX,
    CURVE_GY,
    BN('1')
  );
  static readonly ZERO: JacobianPoint = new JacobianPoint(
    BN('0'),
    BN('1'),
    BN('0')
  );
  static fromAffine(p: Point): JacobianPoint {
    if (!(p instanceof Point)) {
      throw new TypeError('JacobianPoint#fromAffine: expected Point');
    }
    // fromAffine(x:0, y:0) would produce (x:0, y:0, z:1), but we need (x:0, y:1, z:0)
    if (p.equals(Point.ZERO)) return JacobianPoint.ZERO;
    return new JacobianPoint(p.x, p.y, BN('1'));
  }

  /**
   * Compare one point to another.
   */
  equals(other: JacobianPoint): boolean {
    if (!(other instanceof JacobianPoint))
      throw new TypeError('JacobianPoint expected');
    const X1 = this.x;
    const Y1 = this.y;
    const Z1 = this.z;
    const X2 = other.x;
    const Y2 = other.y;
    const Z2 = other.z;

    const Z1Z1 = mod(Z1.mul(Z1));
    const Z2Z2 = mod(Z2.mul(Z2));
    const U1 = mod(X1.mul(Z2Z2));
    const U2 = mod(X2.mul(Z1Z1));
    const S1 = mod(mod(Y1.mul(Z2)).mul(Z2Z2));
    const S2 = mod(mod(Y2.mul(Z1)).mul(Z1Z1));
    return U1.eq(U2) && S1.eq(S2);
  }

  /**
   * Flips point to one corresponding to (x, -y) in Affine coordinates.
   */
  negate(): JacobianPoint {
    return new JacobianPoint(this.x, mod(this.y.opposite()), this.z);
  }

  // Fast algo for doubling 2 Jacobian Points.
  // From: http://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html#doubling-dbl-2007-bl
  // Cost: 1M + 8S + 1*a + 10add + 2*2 + 1*3 + 1*8..
  double(): JacobianPoint {
    const X1 = this.x;
    const Y1 = this.y;
    const Z1 = this.z;

    const XX = mod(X1.mul(X1)); // XX = X1^2
    const YY = mod(Y1.mul(Y1)); // YY = Y1^2
    const YYYY = mod(YY.mul(YY)); // YYYY = YY^2
    const ZZ = mod(Z1.mul(Z1)); // ZZ = Z1^2
    const tmp1 = mod(X1.add(YY)); // (X1+YY)
    const tmp2 = mod(tmp1.mul(tmp1)); // (X1+YY)^2
    const S = mod(BN('2').mul(tmp2.sub(XX).sub(YYYY))); // 2*((X1+YY)^2-XX-YYYY)
    const ZZZZ = mod(ZZ.mul(ZZ)); // ZZ^2
    const M = mod(BN('3').mul(XX).add(CURVE_A.mul(ZZZZ))); // 3*XX+a*ZZ^2
    const MM = mod(M.mul(M)); // M^2
    const T = mod(MM.sub(BN('2').mul(S))); // M^2-2*S
    const X3 = T;
    const Y3 = mod(M.mul(S.sub(T)).sub(BN('8').mul(YYYY))); // M*(S-T)-8*YYYY
    const Y1Z1 = mod(Y1.add(Z1)); // (Y1+Z1)
    const tmp3 = mod(Y1Z1.mul(Y1Z1)); // (Y1+Z1)^2
    const Z3 = mod(tmp3.sub(YY).sub(ZZ)); // (Y1+Z1)^2-YY-ZZ
    return new JacobianPoint(X3, Y3, Z3);
  }

  // Fast algo for adding 2 Jacobian Points.
  // http://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html#addition-add-1998-cmo-2
  // Cost: 12M + 4S + 6add + 1*2
  // Note: 2007 Bernstein-Lange (11M + 5S + 9add + 4*2) is actually 10% slower.
  add(other: JacobianPoint): JacobianPoint {
    if (this.equals(JacobianPoint.ZERO)) return other;
    if (!(other instanceof JacobianPoint))
      throw new TypeError('JacobianPoint expected');
    const X1 = this.x;
    const Y1 = this.y;
    const Z1 = this.z;
    const X2 = other.x;
    const Y2 = other.y;
    const Z2 = other.z;
    if (X2.eq(BN('0')) || Y2.eq(BN('0'))) return this;
    if (X1.eq(BN('0')) || Y1.eq(BN('0'))) return other;
    // We're using same code in equals()
    const Z1Z1 = mod(Z1.mul(Z1)); // Z1Z1 = Z1^2
    const Z2Z2 = mod(Z2.mul(Z2)); // Z2Z2 = Z2^2;
    const U1 = mod(X1.mul(Z2Z2)); // X1 * Z2Z2
    const U2 = mod(X2.mul(Z1Z1)); // X2 * Z1Z1
    const S1 = mod(mod(Y1.mul(Z2)).mul(Z2Z2)); // Y1 * Z2 * Z2Z2
    const S2 = mod(mod(Y2.mul(Z1)).mul(Z1Z1)); // Y2 * Z1 * Z1Z1
    const H = mod(U2.sub(U1)); // H = U2 - U1
    const r = mod(S2.sub(S1)); // S2 - S1
    // H = 0 meaning it's the same point.
    if (H.eq(BN('0'))) {
      if (r.eq(BN('0'))) {
        return this.double();
      } else {
        return JacobianPoint.ZERO;
      }
    }
    const HH = mod(H.mul(H)); // HH = H2
    const HHH = mod(H.mul(HH)); // HHH = H * HH
    const V = mod(U1.mul(HH)); // V = U1 * HH
    const X3 = mod(r.mul(r).sub(HHH).sub(BN('2').mul(V))); // X3 = r^2 - HHH - 2 * V;
    const Y3 = mod(r.mul(V.sub(X3)).sub(S1.mul(HHH))); // Y3 = r * (V - X3) - S1 * HHH;
    const Z3 = mod(Z1.mul(Z2).mul(H)); // Z3 = Z1 * Z2 * H;
    return new JacobianPoint(X3, Y3, Z3);
  }

  subtract(other: JacobianPoint): JacobianPoint {
    return this.add(other.negate());
  }

  // Converts Jacobian point to affine (x, y) coordinates.
  // Can accept precomputed Z^-1 - for example, from invertBatch.
  // (x, y, z) ∋ (x=x/z², y=y/z³)
  // http://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html#scaling-z
  toAffine(): Point {
    const is0 = this.equals(JacobianPoint.ZERO);
    const invZ = is0 ? BN('8') : invert(this.z, CURVE_P); // 8 was chosen arbitrarily
    const iz1 = invZ; // A
    const iz2 = mod(iz1.mul(iz1)); // AA = A^2
    const iz3 = mod(iz2.mul(iz1)); // AAA = A^2 * A = A^3
    const ax = mod(this.x.mul(iz2)); // X3 = X1*AA
    const ay = mod(this.y.mul(iz3)); // Y3 = Y1*AA*A
    const zz = mod(this.z.mul(iz1));
    if (is0) return Point.ZERO;
    if (!zz.eq(BN('1'))) throw new Error('invZ was invalid');
    return new Point(ax, ay);
  }
}

/**
 * Default Point works in default aka affine coordinates: (x, y)
 */
class Point {
  /**
   * Base point aka generator. public_key = Point.BASE * private_key
   */
  static BASE: Point = new Point(CURVE_GX, CURVE_GY);
  /**
   * Identity point aka point at infinity. point = point + zero_point
   */
  static ZERO: Point = new Point(BN('0'), BN('0'));

  constructor(readonly x: BigInt, readonly y: BigInt) {}

  toHex(): string {
    return `04${numTo32bStr(this.x)}${numTo32bStr(this.y)}`;
  }

  toHexX(): string {
    return numTo32bStr(this.x);
  }

  toRawX(): Uint8Array {
    return hexToBytes(this.toHexX());
  }

  equals(other: Point): boolean {
    return this.x.eq(other.x) && this.y.eq(other.y);
  }

  // Returns the same point with inverted `y`
  negate(): Point {
    return new Point(this.x, mod(this.y.opposite()));
  }

  // Adds point to itself
  double(): Point {
    return JacobianPoint.fromAffine(this).double().toAffine();
  }

  // Adds point to other point
  add(other: Point): Point {
    return JacobianPoint.fromAffine(this)
      .add(JacobianPoint.fromAffine(other))
      .toAffine();
  }

  // Subtracts other point from the point
  subtract(other: Point): Point {
    return this.add(other.negate());
  }
}

// Convert between types
// ---------------------

const POW_2_256 = BN(
  '0x10000000000000000000000000000000000000000000000000000000000000000'
);
function numTo32bStr(num: BigInt): string {
  if (!(num instanceof BigInt)) throw new Error('Expected BigInt');
  if (!(BN('0').lte(num) && num.lt(POW_2_256)))
    throw new Error('Expected number < 2^256');
  return num.toString(16).padStart(64, '0');
}

function strip0x(hex: string): string {
  return hex.replace('0x', '');
}

// Caching slows it down 2-3x
function hexToBytes(hex: string): Uint8Array {
  // Stakware has eth-like hexes
  hex = strip0x(hex);
  if (hex.length & 1) hex = `0${hex}`; // padding
  if (typeof hex !== 'string') {
    throw new TypeError('hexToBytes: expected string, got ' + typeof hex);
  }
  if (hex.length % 2)
    throw new Error(`hexToBytes: received invalid unpadded hex ${hex.length}`);
  const array = new Uint8Array(hex.length / 2);
  for (let i = 0; i < array.length; i++) {
    const j = i * 2;
    const hexByte = hex.slice(j, j + 2);
    const byte = i32.parse(hexByte, 16);
    if (Number.isNaN(byte) || byte < 0)
      throw new Error('Invalid byte sequence');
    array[i] = byte;
  }
  return array;
}

const hexes = new Array<number>(256).map((_: number, i: i32) =>
  i.toString(16).padStart(2, '0')
);
function bytesToHex(uint8a: Uint8Array): string {
  if (!(uint8a instanceof Uint8Array)) throw new Error('Expected Uint8Array');
  // pre-caching improves the speed 6x
  let hex = '';
  for (let i = 0; i < uint8a.length; i++) {
    hex += hexes[uint8a[i]];
  }
  return hex;
}

// Regex is not supported
const stripLeadingZeros = (s: string): string => {
  return BigInt.fromString(s, 16).toString(16);
};
const bytesToHexEth = (uint8a: Uint8Array): string =>
  `0x${stripLeadingZeros(bytesToHex(uint8a))}`;

// -------------------------

// Calculates a modulo b
function mod(a: BigInt, b: BigInt = CURVE_P): BigInt {
  const result = a.mod(b);
  return result.gte(BN('0')) ? result : b.add(result);
}

// Inverses number over modulo
function invert(number: BigInt, modulo: BigInt): BigInt {
  if (number.eq(BN('0')) || modulo.lte(BN('0'))) {
    throw new Error(
      `invert: expected positive integers, got n=${number} mod=${modulo}`
    );
  }
  // Eucledian GCD https://brilliant.org/wiki/extended-euclidean-algorithm/
  let a = mod(number, modulo);
  let b = modulo;
  // prettier-ignore
  let x = BN('0'), y = BN('1'), u = BN('1'), v = BN('0');
  while (!a.eq(BN('0'))) {
    const q = b.div(a);
    const r = b.mod(a);
    const m = x.sub(u.mul(q));
    const n = y.sub(v.mul(q));
    (b = a), (a = r), (x = u), (y = v), (u = m), (v = n);
  }
  const gcd = b;
  if (!gcd.eq(BN('1'))) throw new Error('invert: does not exist');
  return mod(x, modulo);
}

// https://docs.starkware.co/starkex/pedersen-hash-function.html
const PEDERSEN_POINTS = [
  new Point(
    BN(
      '2089986280348253421170679821480865132823066470938446095505822317253594081284'
    ),
    BN(
      '1713931329540660377023406109199410414810705867260802078187082345529207694986'
    )
  ),
  new Point(
    BN(
      '996781205833008774514500082376783249102396023663454813447423147977397232763'
    ),
    BN(
      '1668503676786377725805489344771023921079126552019160156920634619255970485781'
    )
  ),
  new Point(
    BN(
      '2251563274489750535117886426533222435294046428347329203627021249169616184184'
    ),
    BN(
      '1798716007562728905295480679789526322175868328062420237419143593021674992973'
    )
  ),
  new Point(
    BN(
      '2138414695194151160943305727036575959195309218611738193261179310511854807447'
    ),
    BN(
      '113410276730064486255102093846540133784865286929052426931474106396135072156'
    )
  ),
  new Point(
    BN(
      '2379962749567351885752724891227938183011949129833673362440656643086021394946'
    ),
    BN(
      '776496453633298175483985398648758586525933812536653089401905292063708816422'
    )
  ),
];
const PEDERSEN_POINTS_JACOBIAN = PEDERSEN_POINTS.map((p: Point) =>
  JacobianPoint.fromAffine(p)
);

function pedersenPrecompute(
  p1: JacobianPoint,
  p2: JacobianPoint
): JacobianPoint[] {
  const out: JacobianPoint[] = [];
  let p = p1;
  for (let i = 0; i < 248; i++) {
    out.push(p);
    p = p.double();
  }
  p = p2;
  for (let i = 0; i < 4; i++) {
    out.push(p);
    p = p.double();
  }
  return out;
}
const PEDERSEN_POINTS1: JacobianPoint[] = pedersenPrecompute(
  PEDERSEN_POINTS_JACOBIAN[1],
  PEDERSEN_POINTS_JACOBIAN[2]
);

const PEDERSEN_POINTS2: JacobianPoint[] = pedersenPrecompute(
  PEDERSEN_POINTS_JACOBIAN[3],
  PEDERSEN_POINTS_JACOBIAN[4]
);

function pedersenArg(value: BigInt): BigInt {
  // [0..Fp)
  if (BN('0').gt(value) || value.gte(CURVE_P)) {
    throw new Error(`value should be 0<=ARG<CURVE.P: ${value}`);
  }
  return value;
}

function pedersenSingle(
  point: JacobianPoint,
  value: BigInt,
  constants: JacobianPoint[]
): JacobianPoint {
  let x = pedersenArg(value);

  for (let j = 0; j < 252; j++) {
    const pt = constants[j];
    if (pt.x.eq(point.x)) throw new Error('Same point');
    // Hacky fix for bigint comparision issue
    if (BigInt.from(x.bitwiseAnd(BN('1')).toString()).ne(BN('0'))) {
      point = point.add(pt);
    }
    x = x.rightShift(1);
  }
  return point;
}

// shift_point + x_low * P_0 + x_high * P1 + y_low * P2  + y_high * P3
export function pedersen(x: string, y: string): string {
  let point = PEDERSEN_POINTS_JACOBIAN[0];
  point = pedersenSingle(point, BN(x), PEDERSEN_POINTS1);
  point = pedersenSingle(point, BN(y), PEDERSEN_POINTS2);

  return bytesToHexEth(point.toAffine().toRawX());
}

export function computeHashOnElements(data: string[]): string {
  return data.concat([data.length.toString()]).reduce((x, y) => {
    return pedersen(x, y);
  }, '0');
}
