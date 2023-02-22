import {
  PEDERSEN_TEST_VALUES,
  COMPUTE_HASH_ON_ELEMENTS_TEST_VALUES,
} from './utils';
import { pedersen, computeHashOnElements } from '../index';

describe('Pedersen Hashing', () => {
  it('should compute the correct pedersen hash', () => {
    for (let i = 0; i < PEDERSEN_TEST_VALUES.length; i += 3) {
      const x = PEDERSEN_TEST_VALUES[i];
      const y = PEDERSEN_TEST_VALUES[i + 1];
      const output = PEDERSEN_TEST_VALUES[i + 2];
      expect<String>(pedersen(x, y)).toStrictEqual(output);
    }
  });

  it('should compute the correct pedersen hash with many elements', () => {
    for (let i = 0; i < COMPUTE_HASH_ON_ELEMENTS_TEST_VALUES.length; i += 4) {
      const a = COMPUTE_HASH_ON_ELEMENTS_TEST_VALUES[i];
      const b = COMPUTE_HASH_ON_ELEMENTS_TEST_VALUES[i + 1];
      const c = COMPUTE_HASH_ON_ELEMENTS_TEST_VALUES[i + 2];
      const output = COMPUTE_HASH_ON_ELEMENTS_TEST_VALUES[i + 3];
      expect<String>(computeHashOnElements([a, b, c])).toStrictEqual(output);
    }
  });
});
