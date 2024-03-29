/******************************************************************************
Copyright (c) Microsoft Corporation.

Permission to use, copy, modify, and/or distribute this software for any
purpose with or without fee is hereby granted.

THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH
REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY
AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT,
INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM
LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
PERFORMANCE OF THIS SOFTWARE.
***************************************************************************** */
/* global Reflect, Promise, SuppressedError, Symbol */


function __awaiter(thisArg, _arguments, P, generator) {
    function adopt(value) { return value instanceof P ? value : new P(function (resolve) { resolve(value); }); }
    return new (P || (P = Promise))(function (resolve, reject) {
        function fulfilled(value) { try { step(generator.next(value)); } catch (e) { reject(e); } }
        function rejected(value) { try { step(generator["throw"](value)); } catch (e) { reject(e); } }
        function step(result) { result.done ? resolve(result.value) : adopt(result.value).then(fulfilled, rejected); }
        step((generator = generator.apply(thisArg, _arguments || [])).next());
    });
}

typeof SuppressedError === "function" ? SuppressedError : function (error, suppressed, message) {
    var e = new Error(message);
    return e.name = "SuppressedError", e.error = error, e.suppressed = suppressed, e;
};

const converterToEng = {
    кс: "x",
    щ: "sch",
    ч: "ch",
    ж: "zh",
    х: "kh",
    ш: "sh",
    ю: "yu",
    я: "ya",
    а: "a",
    б: "b",
    в: "v",
    г: "g",
    д: "d",
    е: "e",
    ё: "e",
    з: "z",
    и: "i",
    й: "y",
    к: "k",
    л: "l",
    м: "m",
    н: "n",
    о: "o",
    п: "p",
    р: "r",
    с: "s",
    т: "t",
    у: "u",
    ф: "f",
    ц: "c",
    ь: "",
    ы: "y",
    ъ: "",
    э: "e",
};
function getOtherLanguageStr(str, converter) {
    Object.keys(converter).map((sub) => {
        str = str.replace(new RegExp(sub, "g"), converter[sub]);
    });
    return str;
}

function getEnglishStr(str) {
    return getOtherLanguageStr(str, converterToEng);
}

const converterToRus = {};
Object.entries(converterToEng).forEach(([rus, eng]) => {
    if (!(converterToRus[eng] || eng === ""))
        converterToRus[eng] = rus;
});
function getRussianStr(str) {
    return getOtherLanguageStr(str, converterToRus);
}

class StringRewriter {
    constructor() {
        this._keyboardRuToEn = new Map([
            ["й", "q"],
            ["ц", "w"],
            ["у", "e"],
            ["к", "r"],
            ["е", "t"],
            ["н", "y"],
            ["г", "u"],
            ["ш", "i"],
            ["щ", "o"],
            ["з", "p"],
            ["х", "["],
            ["ъ", "]"],
            ["ф", "a"],
            ["ы", "s"],
            ["в", "d"],
            ["а", "f"],
            ["п", "g"],
            ["р", "h"],
            ["о", "j"],
            ["л", "k"],
            ["д", "l"],
            ["ж", ";"],
            ["э", "'"],
            ["я", "z"],
            ["ч", "x"],
            ["с", "c"],
            ["м", "v"],
            ["и", "b"],
            ["т", "n"],
            ["ь", "m"],
            ["б", ","],
            ["ю", "."],
        ]);
    }
    changeTextLayout(stringToFix) {
        let wrongLayoutQuery = "";
        stringToFix.toLowerCase();
        for (let i = 0; i < stringToFix.length; i++) {
            const letter = stringToFix[i];
            const char = this._keyboardRuToEn.get(letter) || this._getKeyByValue(letter) || letter;
            wrongLayoutQuery += char;
        }
        return wrongLayoutQuery;
    }
    changeRussianToEnglishTransliteration(stringToFix) {
        return getEnglishStr(stringToFix);
    }
    changeEnglishToRussianTransliteration(stringToFix) {
        return getRussianStr(stringToFix);
    }
    _getKeyByValue(letter) {
        const keys = this._keyboardRuToEn.keys();
        for (const key of keys)
            if (this._keyboardRuToEn.get(key) === letter)
                return key;
        return null;
    }
}

const multiLayoutSearcher = (searcher) => {
    return (query) => __awaiter(void 0, void 0, void 0, function* () {
        let result = yield searcher(query);
        if (result)
            return result;
        const stringRewriter = new StringRewriter();
        const wrongLayoutQuery = stringRewriter.changeTextLayout(query);
        result = yield searcher(wrongLayoutQuery);
        if (result)
            return result;
        const RuToEnRev = stringRewriter.changeRussianToEnglishTransliteration(query);
        result = yield searcher(RuToEnRev);
        if (result)
            return result;
        const EnToRuRev = stringRewriter.changeEnglishToRussianTransliteration(query);
        result = yield searcher(EnToRuRev);
        return result;
    });
};

var ArticleType;
(function (ArticleType) {
    ArticleType["article"] = "article";
    ArticleType["category"] = "category";
})(ArticleType || (ArticleType = {}));

var FileStatus;
(function (FileStatus) {
    FileStatus["modified"] = "modified";
    FileStatus["delete"] = "delete";
    FileStatus["new"] = "new";
    FileStatus["rename"] = "rename";
    FileStatus["conflict"] = "conflict";
    FileStatus["current"] = "current";
})(FileStatus || (FileStatus = {}));

class Plugin {
    constructor(_app) {
        this._app = _app;
        this._commandConfigs = [];
    }
    get commandConfigs() {
        return this._commandConfigs;
    }
    addCommand(command) {
        this._commandConfigs.push(command);
    }
}

/**
 * Fuse.js v7.0.0 - Lightweight fuzzy-search (http://fusejs.io)
 *
 * Copyright (c) 2023 Kiro Risk (http://kiro.me)
 * All Rights Reserved. Apache Software License 2.0
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 */

function isArray(value) {
  return !Array.isArray
    ? getTag(value) === '[object Array]'
    : Array.isArray(value)
}

// Adapted from: https://github.com/lodash/lodash/blob/master/.internal/baseToString.js
const INFINITY = 1 / 0;
function baseToString(value) {
  // Exit early for strings to avoid a performance hit in some environments.
  if (typeof value == 'string') {
    return value
  }
  let result = value + '';
  return result == '0' && 1 / value == -INFINITY ? '-0' : result
}

function toString(value) {
  return value == null ? '' : baseToString(value)
}

function isString(value) {
  return typeof value === 'string'
}

function isNumber(value) {
  return typeof value === 'number'
}

// Adapted from: https://github.com/lodash/lodash/blob/master/isBoolean.js
function isBoolean(value) {
  return (
    value === true ||
    value === false ||
    (isObjectLike(value) && getTag(value) == '[object Boolean]')
  )
}

function isObject(value) {
  return typeof value === 'object'
}

// Checks if `value` is object-like.
function isObjectLike(value) {
  return isObject(value) && value !== null
}

function isDefined(value) {
  return value !== undefined && value !== null
}

function isBlank(value) {
  return !value.trim().length
}

// Gets the `toStringTag` of `value`.
// Adapted from: https://github.com/lodash/lodash/blob/master/.internal/getTag.js
function getTag(value) {
  return value == null
    ? value === undefined
      ? '[object Undefined]'
      : '[object Null]'
    : Object.prototype.toString.call(value)
}

const EXTENDED_SEARCH_UNAVAILABLE = 'Extended search is not available';

const INCORRECT_INDEX_TYPE = "Incorrect 'index' type";

const LOGICAL_SEARCH_INVALID_QUERY_FOR_KEY = (key) =>
  `Invalid value for key ${key}`;

const PATTERN_LENGTH_TOO_LARGE = (max) =>
  `Pattern length exceeds max of ${max}.`;

const MISSING_KEY_PROPERTY = (name) => `Missing ${name} property in key`;

const INVALID_KEY_WEIGHT_VALUE = (key) =>
  `Property 'weight' in key '${key}' must be a positive integer`;

const hasOwn = Object.prototype.hasOwnProperty;

class KeyStore {
  constructor(keys) {
    this._keys = [];
    this._keyMap = {};

    let totalWeight = 0;

    keys.forEach((key) => {
      let obj = createKey(key);

      this._keys.push(obj);
      this._keyMap[obj.id] = obj;

      totalWeight += obj.weight;
    });

    // Normalize weights so that their sum is equal to 1
    this._keys.forEach((key) => {
      key.weight /= totalWeight;
    });
  }
  get(keyId) {
    return this._keyMap[keyId]
  }
  keys() {
    return this._keys
  }
  toJSON() {
    return JSON.stringify(this._keys)
  }
}

function createKey(key) {
  let path = null;
  let id = null;
  let src = null;
  let weight = 1;
  let getFn = null;

  if (isString(key) || isArray(key)) {
    src = key;
    path = createKeyPath(key);
    id = createKeyId(key);
  } else {
    if (!hasOwn.call(key, 'name')) {
      throw new Error(MISSING_KEY_PROPERTY('name'))
    }

    const name = key.name;
    src = name;

    if (hasOwn.call(key, 'weight')) {
      weight = key.weight;

      if (weight <= 0) {
        throw new Error(INVALID_KEY_WEIGHT_VALUE(name))
      }
    }

    path = createKeyPath(name);
    id = createKeyId(name);
    getFn = key.getFn;
  }

  return { path, id, weight, src, getFn }
}

function createKeyPath(key) {
  return isArray(key) ? key : key.split('.')
}

function createKeyId(key) {
  return isArray(key) ? key.join('.') : key
}

function get$3(obj, path) {
  let list = [];
  let arr = false;

  const deepGet = (obj, path, index) => {
    if (!isDefined(obj)) {
      return
    }
    if (!path[index]) {
      // If there's no path left, we've arrived at the object we care about.
      list.push(obj);
    } else {
      let key = path[index];

      const value = obj[key];

      if (!isDefined(value)) {
        return
      }

      // If we're at the last value in the path, and if it's a string/number/bool,
      // add it to the list
      if (
        index === path.length - 1 &&
        (isString(value) || isNumber(value) || isBoolean(value))
      ) {
        list.push(toString(value));
      } else if (isArray(value)) {
        arr = true;
        // Search each item in the array.
        for (let i = 0, len = value.length; i < len; i += 1) {
          deepGet(value[i], path, index + 1);
        }
      } else if (path.length) {
        // An object. Recurse further.
        deepGet(value, path, index + 1);
      }
    }
  };

  // Backwards compatibility (since path used to be a string)
  deepGet(obj, isString(path) ? path.split('.') : path, 0);

  return arr ? list : list[0]
}

const MatchOptions = {
  // Whether the matches should be included in the result set. When `true`, each record in the result
  // set will include the indices of the matched characters.
  // These can consequently be used for highlighting purposes.
  includeMatches: false,
  // When `true`, the matching function will continue to the end of a search pattern even if
  // a perfect match has already been located in the string.
  findAllMatches: false,
  // Minimum number of characters that must be matched before a result is considered a match
  minMatchCharLength: 1
};

const BasicOptions = {
  // When `true`, the algorithm continues searching to the end of the input even if a perfect
  // match is found before the end of the same input.
  isCaseSensitive: false,
  // When true, the matching function will continue to the end of a search pattern even if
  includeScore: false,
  // List of properties that will be searched. This also supports nested properties.
  keys: [],
  // Whether to sort the result list, by score
  shouldSort: true,
  // Default sort function: sort by ascending score, ascending index
  sortFn: (a, b) =>
    a.score === b.score ? (a.idx < b.idx ? -1 : 1) : a.score < b.score ? -1 : 1
};

const FuzzyOptions = {
  // Approximately where in the text is the pattern expected to be found?
  location: 0,
  // At what point does the match algorithm give up. A threshold of '0.0' requires a perfect match
  // (of both letters and location), a threshold of '1.0' would match anything.
  threshold: 0.6,
  // Determines how close the match must be to the fuzzy location (specified above).
  // An exact letter match which is 'distance' characters away from the fuzzy location
  // would score as a complete mismatch. A distance of '0' requires the match be at
  // the exact location specified, a threshold of '1000' would require a perfect match
  // to be within 800 characters of the fuzzy location to be found using a 0.8 threshold.
  distance: 100
};

const AdvancedOptions = {
  // When `true`, it enables the use of unix-like search commands
  useExtendedSearch: false,
  // The get function to use when fetching an object's properties.
  // The default will search nested paths *ie foo.bar.baz*
  getFn: get$3,
  // When `true`, search will ignore `location` and `distance`, so it won't matter
  // where in the string the pattern appears.
  // More info: https://fusejs.io/concepts/scoring-theory.html#fuzziness-score
  ignoreLocation: false,
  // When `true`, the calculation for the relevance score (used for sorting) will
  // ignore the field-length norm.
  // More info: https://fusejs.io/concepts/scoring-theory.html#field-length-norm
  ignoreFieldNorm: false,
  // The weight to determine how much field length norm effects scoring.
  fieldNormWeight: 1
};

var Config = {
  ...BasicOptions,
  ...MatchOptions,
  ...FuzzyOptions,
  ...AdvancedOptions
};

const SPACE = /[^ ]+/g;

// Field-length norm: the shorter the field, the higher the weight.
// Set to 3 decimals to reduce index size.
function norm(weight = 1, mantissa = 3) {
  const cache = new Map();
  const m = Math.pow(10, mantissa);

  return {
    get(value) {
      const numTokens = value.match(SPACE).length;

      if (cache.has(numTokens)) {
        return cache.get(numTokens)
      }

      // Default function is 1/sqrt(x), weight makes that variable
      const norm = 1 / Math.pow(numTokens, 0.5 * weight);

      // In place of `toFixed(mantissa)`, for faster computation
      const n = parseFloat(Math.round(norm * m) / m);

      cache.set(numTokens, n);

      return n
    },
    clear() {
      cache.clear();
    }
  }
}

class FuseIndex {
  constructor({
    getFn = Config.getFn,
    fieldNormWeight = Config.fieldNormWeight
  } = {}) {
    this.norm = norm(fieldNormWeight, 3);
    this.getFn = getFn;
    this.isCreated = false;

    this.setIndexRecords();
  }
  setSources(docs = []) {
    this.docs = docs;
  }
  setIndexRecords(records = []) {
    this.records = records;
  }
  setKeys(keys = []) {
    this.keys = keys;
    this._keysMap = {};
    keys.forEach((key, idx) => {
      this._keysMap[key.id] = idx;
    });
  }
  create() {
    if (this.isCreated || !this.docs.length) {
      return
    }

    this.isCreated = true;

    // List is Array<String>
    if (isString(this.docs[0])) {
      this.docs.forEach((doc, docIndex) => {
        this._addString(doc, docIndex);
      });
    } else {
      // List is Array<Object>
      this.docs.forEach((doc, docIndex) => {
        this._addObject(doc, docIndex);
      });
    }

    this.norm.clear();
  }
  // Adds a doc to the end of the index
  add(doc) {
    const idx = this.size();

    if (isString(doc)) {
      this._addString(doc, idx);
    } else {
      this._addObject(doc, idx);
    }
  }
  // Removes the doc at the specified index of the index
  removeAt(idx) {
    this.records.splice(idx, 1);

    // Change ref index of every subsquent doc
    for (let i = idx, len = this.size(); i < len; i += 1) {
      this.records[i].i -= 1;
    }
  }
  getValueForItemAtKeyId(item, keyId) {
    return item[this._keysMap[keyId]]
  }
  size() {
    return this.records.length
  }
  _addString(doc, docIndex) {
    if (!isDefined(doc) || isBlank(doc)) {
      return
    }

    let record = {
      v: doc,
      i: docIndex,
      n: this.norm.get(doc)
    };

    this.records.push(record);
  }
  _addObject(doc, docIndex) {
    let record = { i: docIndex, $: {} };

    // Iterate over every key (i.e, path), and fetch the value at that key
    this.keys.forEach((key, keyIndex) => {
      let value = key.getFn ? key.getFn(doc) : this.getFn(doc, key.path);

      if (!isDefined(value)) {
        return
      }

      if (isArray(value)) {
        let subRecords = [];
        const stack = [{ nestedArrIndex: -1, value }];

        while (stack.length) {
          const { nestedArrIndex, value } = stack.pop();

          if (!isDefined(value)) {
            continue
          }

          if (isString(value) && !isBlank(value)) {
            let subRecord = {
              v: value,
              i: nestedArrIndex,
              n: this.norm.get(value)
            };

            subRecords.push(subRecord);
          } else if (isArray(value)) {
            value.forEach((item, k) => {
              stack.push({
                nestedArrIndex: k,
                value: item
              });
            });
          } else ;
        }
        record.$[keyIndex] = subRecords;
      } else if (isString(value) && !isBlank(value)) {
        let subRecord = {
          v: value,
          n: this.norm.get(value)
        };

        record.$[keyIndex] = subRecord;
      }
    });

    this.records.push(record);
  }
  toJSON() {
    return {
      keys: this.keys,
      records: this.records
    }
  }
}

function createIndex(
  keys,
  docs,
  { getFn = Config.getFn, fieldNormWeight = Config.fieldNormWeight } = {}
) {
  const myIndex = new FuseIndex({ getFn, fieldNormWeight });
  myIndex.setKeys(keys.map(createKey));
  myIndex.setSources(docs);
  myIndex.create();
  return myIndex
}

function parseIndex(
  data,
  { getFn = Config.getFn, fieldNormWeight = Config.fieldNormWeight } = {}
) {
  const { keys, records } = data;
  const myIndex = new FuseIndex({ getFn, fieldNormWeight });
  myIndex.setKeys(keys);
  myIndex.setIndexRecords(records);
  return myIndex
}

function computeScore$1(
  pattern,
  {
    errors = 0,
    currentLocation = 0,
    expectedLocation = 0,
    distance = Config.distance,
    ignoreLocation = Config.ignoreLocation
  } = {}
) {
  const accuracy = errors / pattern.length;

  if (ignoreLocation) {
    return accuracy
  }

  const proximity = Math.abs(expectedLocation - currentLocation);

  if (!distance) {
    // Dodge divide by zero error.
    return proximity ? 1.0 : accuracy
  }

  return accuracy + proximity / distance
}

function convertMaskToIndices(
  matchmask = [],
  minMatchCharLength = Config.minMatchCharLength
) {
  let indices = [];
  let start = -1;
  let end = -1;
  let i = 0;

  for (let len = matchmask.length; i < len; i += 1) {
    let match = matchmask[i];
    if (match && start === -1) {
      start = i;
    } else if (!match && start !== -1) {
      end = i - 1;
      if (end - start + 1 >= minMatchCharLength) {
        indices.push([start, end]);
      }
      start = -1;
    }
  }

  // (i-1 - start) + 1 => i - start
  if (matchmask[i - 1] && i - start >= minMatchCharLength) {
    indices.push([start, i - 1]);
  }

  return indices
}

// Machine word size
const MAX_BITS = 32;

function search(
  text,
  pattern,
  patternAlphabet,
  {
    location = Config.location,
    distance = Config.distance,
    threshold = Config.threshold,
    findAllMatches = Config.findAllMatches,
    minMatchCharLength = Config.minMatchCharLength,
    includeMatches = Config.includeMatches,
    ignoreLocation = Config.ignoreLocation
  } = {}
) {
  if (pattern.length > MAX_BITS) {
    throw new Error(PATTERN_LENGTH_TOO_LARGE(MAX_BITS))
  }

  const patternLen = pattern.length;
  // Set starting location at beginning text and initialize the alphabet.
  const textLen = text.length;
  // Handle the case when location > text.length
  const expectedLocation = Math.max(0, Math.min(location, textLen));
  // Highest score beyond which we give up.
  let currentThreshold = threshold;
  // Is there a nearby exact match? (speedup)
  let bestLocation = expectedLocation;

  // Performance: only computer matches when the minMatchCharLength > 1
  // OR if `includeMatches` is true.
  const computeMatches = minMatchCharLength > 1 || includeMatches;
  // A mask of the matches, used for building the indices
  const matchMask = computeMatches ? Array(textLen) : [];

  let index;

  // Get all exact matches, here for speed up
  while ((index = text.indexOf(pattern, bestLocation)) > -1) {
    let score = computeScore$1(pattern, {
      currentLocation: index,
      expectedLocation,
      distance,
      ignoreLocation
    });

    currentThreshold = Math.min(score, currentThreshold);
    bestLocation = index + patternLen;

    if (computeMatches) {
      let i = 0;
      while (i < patternLen) {
        matchMask[index + i] = 1;
        i += 1;
      }
    }
  }

  // Reset the best location
  bestLocation = -1;

  let lastBitArr = [];
  let finalScore = 1;
  let binMax = patternLen + textLen;

  const mask = 1 << (patternLen - 1);

  for (let i = 0; i < patternLen; i += 1) {
    // Scan for the best match; each iteration allows for one more error.
    // Run a binary search to determine how far from the match location we can stray
    // at this error level.
    let binMin = 0;
    let binMid = binMax;

    while (binMin < binMid) {
      const score = computeScore$1(pattern, {
        errors: i,
        currentLocation: expectedLocation + binMid,
        expectedLocation,
        distance,
        ignoreLocation
      });

      if (score <= currentThreshold) {
        binMin = binMid;
      } else {
        binMax = binMid;
      }

      binMid = Math.floor((binMax - binMin) / 2 + binMin);
    }

    // Use the result from this iteration as the maximum for the next.
    binMax = binMid;

    let start = Math.max(1, expectedLocation - binMid + 1);
    let finish = findAllMatches
      ? textLen
      : Math.min(expectedLocation + binMid, textLen) + patternLen;

    // Initialize the bit array
    let bitArr = Array(finish + 2);

    bitArr[finish + 1] = (1 << i) - 1;

    for (let j = finish; j >= start; j -= 1) {
      let currentLocation = j - 1;
      let charMatch = patternAlphabet[text.charAt(currentLocation)];

      if (computeMatches) {
        // Speed up: quick bool to int conversion (i.e, `charMatch ? 1 : 0`)
        matchMask[currentLocation] = +!!charMatch;
      }

      // First pass: exact match
      bitArr[j] = ((bitArr[j + 1] << 1) | 1) & charMatch;

      // Subsequent passes: fuzzy match
      if (i) {
        bitArr[j] |=
          ((lastBitArr[j + 1] | lastBitArr[j]) << 1) | 1 | lastBitArr[j + 1];
      }

      if (bitArr[j] & mask) {
        finalScore = computeScore$1(pattern, {
          errors: i,
          currentLocation,
          expectedLocation,
          distance,
          ignoreLocation
        });

        // This match will almost certainly be better than any existing match.
        // But check anyway.
        if (finalScore <= currentThreshold) {
          // Indeed it is
          currentThreshold = finalScore;
          bestLocation = currentLocation;

          // Already passed `loc`, downhill from here on in.
          if (bestLocation <= expectedLocation) {
            break
          }

          // When passing `bestLocation`, don't exceed our current distance from `expectedLocation`.
          start = Math.max(1, 2 * expectedLocation - bestLocation);
        }
      }
    }

    // No hope for a (better) match at greater error levels.
    const score = computeScore$1(pattern, {
      errors: i + 1,
      currentLocation: expectedLocation,
      expectedLocation,
      distance,
      ignoreLocation
    });

    if (score > currentThreshold) {
      break
    }

    lastBitArr = bitArr;
  }

  const result = {
    isMatch: bestLocation >= 0,
    // Count exact matches (those with a score of 0) to be "almost" exact
    score: Math.max(0.001, finalScore)
  };

  if (computeMatches) {
    const indices = convertMaskToIndices(matchMask, minMatchCharLength);
    if (!indices.length) {
      result.isMatch = false;
    } else if (includeMatches) {
      result.indices = indices;
    }
  }

  return result
}

function createPatternAlphabet(pattern) {
  let mask = {};

  for (let i = 0, len = pattern.length; i < len; i += 1) {
    const char = pattern.charAt(i);
    mask[char] = (mask[char] || 0) | (1 << (len - i - 1));
  }

  return mask
}

class BitapSearch {
  constructor(
    pattern,
    {
      location = Config.location,
      threshold = Config.threshold,
      distance = Config.distance,
      includeMatches = Config.includeMatches,
      findAllMatches = Config.findAllMatches,
      minMatchCharLength = Config.minMatchCharLength,
      isCaseSensitive = Config.isCaseSensitive,
      ignoreLocation = Config.ignoreLocation
    } = {}
  ) {
    this.options = {
      location,
      threshold,
      distance,
      includeMatches,
      findAllMatches,
      minMatchCharLength,
      isCaseSensitive,
      ignoreLocation
    };

    this.pattern = isCaseSensitive ? pattern : pattern.toLowerCase();

    this.chunks = [];

    if (!this.pattern.length) {
      return
    }

    const addChunk = (pattern, startIndex) => {
      this.chunks.push({
        pattern,
        alphabet: createPatternAlphabet(pattern),
        startIndex
      });
    };

    const len = this.pattern.length;

    if (len > MAX_BITS) {
      let i = 0;
      const remainder = len % MAX_BITS;
      const end = len - remainder;

      while (i < end) {
        addChunk(this.pattern.substr(i, MAX_BITS), i);
        i += MAX_BITS;
      }

      if (remainder) {
        const startIndex = len - MAX_BITS;
        addChunk(this.pattern.substr(startIndex), startIndex);
      }
    } else {
      addChunk(this.pattern, 0);
    }
  }

  searchIn(text) {
    const { isCaseSensitive, includeMatches } = this.options;

    if (!isCaseSensitive) {
      text = text.toLowerCase();
    }

    // Exact match
    if (this.pattern === text) {
      let result = {
        isMatch: true,
        score: 0
      };

      if (includeMatches) {
        result.indices = [[0, text.length - 1]];
      }

      return result
    }

    // Otherwise, use Bitap algorithm
    const {
      location,
      distance,
      threshold,
      findAllMatches,
      minMatchCharLength,
      ignoreLocation
    } = this.options;

    let allIndices = [];
    let totalScore = 0;
    let hasMatches = false;

    this.chunks.forEach(({ pattern, alphabet, startIndex }) => {
      const { isMatch, score, indices } = search(text, pattern, alphabet, {
        location: location + startIndex,
        distance,
        threshold,
        findAllMatches,
        minMatchCharLength,
        includeMatches,
        ignoreLocation
      });

      if (isMatch) {
        hasMatches = true;
      }

      totalScore += score;

      if (isMatch && indices) {
        allIndices = [...allIndices, ...indices];
      }
    });

    let result = {
      isMatch: hasMatches,
      score: hasMatches ? totalScore / this.chunks.length : 1
    };

    if (hasMatches && includeMatches) {
      result.indices = allIndices;
    }

    return result
  }
}

class BaseMatch {
  constructor(pattern) {
    this.pattern = pattern;
  }
  static isMultiMatch(pattern) {
    return getMatch(pattern, this.multiRegex)
  }
  static isSingleMatch(pattern) {
    return getMatch(pattern, this.singleRegex)
  }
  search(/*text*/) {}
}

function getMatch(pattern, exp) {
  const matches = pattern.match(exp);
  return matches ? matches[1] : null
}

// Token: 'file

class ExactMatch extends BaseMatch {
  constructor(pattern) {
    super(pattern);
  }
  static get type() {
    return 'exact'
  }
  static get multiRegex() {
    return /^="(.*)"$/
  }
  static get singleRegex() {
    return /^=(.*)$/
  }
  search(text) {
    const isMatch = text === this.pattern;

    return {
      isMatch,
      score: isMatch ? 0 : 1,
      indices: [0, this.pattern.length - 1]
    }
  }
}

// Token: !fire

class InverseExactMatch extends BaseMatch {
  constructor(pattern) {
    super(pattern);
  }
  static get type() {
    return 'inverse-exact'
  }
  static get multiRegex() {
    return /^!"(.*)"$/
  }
  static get singleRegex() {
    return /^!(.*)$/
  }
  search(text) {
    const index = text.indexOf(this.pattern);
    const isMatch = index === -1;

    return {
      isMatch,
      score: isMatch ? 0 : 1,
      indices: [0, text.length - 1]
    }
  }
}

// Token: ^file

class PrefixExactMatch extends BaseMatch {
  constructor(pattern) {
    super(pattern);
  }
  static get type() {
    return 'prefix-exact'
  }
  static get multiRegex() {
    return /^\^"(.*)"$/
  }
  static get singleRegex() {
    return /^\^(.*)$/
  }
  search(text) {
    const isMatch = text.startsWith(this.pattern);

    return {
      isMatch,
      score: isMatch ? 0 : 1,
      indices: [0, this.pattern.length - 1]
    }
  }
}

// Token: !^fire

class InversePrefixExactMatch extends BaseMatch {
  constructor(pattern) {
    super(pattern);
  }
  static get type() {
    return 'inverse-prefix-exact'
  }
  static get multiRegex() {
    return /^!\^"(.*)"$/
  }
  static get singleRegex() {
    return /^!\^(.*)$/
  }
  search(text) {
    const isMatch = !text.startsWith(this.pattern);

    return {
      isMatch,
      score: isMatch ? 0 : 1,
      indices: [0, text.length - 1]
    }
  }
}

// Token: .file$

class SuffixExactMatch extends BaseMatch {
  constructor(pattern) {
    super(pattern);
  }
  static get type() {
    return 'suffix-exact'
  }
  static get multiRegex() {
    return /^"(.*)"\$$/
  }
  static get singleRegex() {
    return /^(.*)\$$/
  }
  search(text) {
    const isMatch = text.endsWith(this.pattern);

    return {
      isMatch,
      score: isMatch ? 0 : 1,
      indices: [text.length - this.pattern.length, text.length - 1]
    }
  }
}

// Token: !.file$

class InverseSuffixExactMatch extends BaseMatch {
  constructor(pattern) {
    super(pattern);
  }
  static get type() {
    return 'inverse-suffix-exact'
  }
  static get multiRegex() {
    return /^!"(.*)"\$$/
  }
  static get singleRegex() {
    return /^!(.*)\$$/
  }
  search(text) {
    const isMatch = !text.endsWith(this.pattern);
    return {
      isMatch,
      score: isMatch ? 0 : 1,
      indices: [0, text.length - 1]
    }
  }
}

class FuzzyMatch extends BaseMatch {
  constructor(
    pattern,
    {
      location = Config.location,
      threshold = Config.threshold,
      distance = Config.distance,
      includeMatches = Config.includeMatches,
      findAllMatches = Config.findAllMatches,
      minMatchCharLength = Config.minMatchCharLength,
      isCaseSensitive = Config.isCaseSensitive,
      ignoreLocation = Config.ignoreLocation
    } = {}
  ) {
    super(pattern);
    this._bitapSearch = new BitapSearch(pattern, {
      location,
      threshold,
      distance,
      includeMatches,
      findAllMatches,
      minMatchCharLength,
      isCaseSensitive,
      ignoreLocation
    });
  }
  static get type() {
    return 'fuzzy'
  }
  static get multiRegex() {
    return /^"(.*)"$/
  }
  static get singleRegex() {
    return /^(.*)$/
  }
  search(text) {
    return this._bitapSearch.searchIn(text)
  }
}

// Token: 'file

class IncludeMatch extends BaseMatch {
  constructor(pattern) {
    super(pattern);
  }
  static get type() {
    return 'include'
  }
  static get multiRegex() {
    return /^'"(.*)"$/
  }
  static get singleRegex() {
    return /^'(.*)$/
  }
  search(text) {
    let location = 0;
    let index;

    const indices = [];
    const patternLen = this.pattern.length;

    // Get all exact matches
    while ((index = text.indexOf(this.pattern, location)) > -1) {
      location = index + patternLen;
      indices.push([index, location - 1]);
    }

    const isMatch = !!indices.length;

    return {
      isMatch,
      score: isMatch ? 0 : 1,
      indices
    }
  }
}

// ❗Order is important. DO NOT CHANGE.
const searchers = [
  ExactMatch,
  IncludeMatch,
  PrefixExactMatch,
  InversePrefixExactMatch,
  InverseSuffixExactMatch,
  SuffixExactMatch,
  InverseExactMatch,
  FuzzyMatch
];

const searchersLen = searchers.length;

// Regex to split by spaces, but keep anything in quotes together
const SPACE_RE = / +(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)/;
const OR_TOKEN = '|';

// Return a 2D array representation of the query, for simpler parsing.
// Example:
// "^core go$ | rb$ | py$ xy$" => [["^core", "go$"], ["rb$"], ["py$", "xy$"]]
function parseQuery(pattern, options = {}) {
  return pattern.split(OR_TOKEN).map((item) => {
    let query = item
      .trim()
      .split(SPACE_RE)
      .filter((item) => item && !!item.trim());

    let results = [];
    for (let i = 0, len = query.length; i < len; i += 1) {
      const queryItem = query[i];

      // 1. Handle multiple query match (i.e, once that are quoted, like `"hello world"`)
      let found = false;
      let idx = -1;
      while (!found && ++idx < searchersLen) {
        const searcher = searchers[idx];
        let token = searcher.isMultiMatch(queryItem);
        if (token) {
          results.push(new searcher(token, options));
          found = true;
        }
      }

      if (found) {
        continue
      }

      // 2. Handle single query matches (i.e, once that are *not* quoted)
      idx = -1;
      while (++idx < searchersLen) {
        const searcher = searchers[idx];
        let token = searcher.isSingleMatch(queryItem);
        if (token) {
          results.push(new searcher(token, options));
          break
        }
      }
    }

    return results
  })
}

// These extended matchers can return an array of matches, as opposed
// to a singl match
const MultiMatchSet = new Set([FuzzyMatch.type, IncludeMatch.type]);

/**
 * Command-like searching
 * ======================
 *
 * Given multiple search terms delimited by spaces.e.g. `^jscript .python$ ruby !java`,
 * search in a given text.
 *
 * Search syntax:
 *
 * | Token       | Match type                 | Description                            |
 * | ----------- | -------------------------- | -------------------------------------- |
 * | `jscript`   | fuzzy-match                | Items that fuzzy match `jscript`       |
 * | `=scheme`   | exact-match                | Items that are `scheme`                |
 * | `'python`   | include-match              | Items that include `python`            |
 * | `!ruby`     | inverse-exact-match        | Items that do not include `ruby`       |
 * | `^java`     | prefix-exact-match         | Items that start with `java`           |
 * | `!^earlang` | inverse-prefix-exact-match | Items that do not start with `earlang` |
 * | `.js$`      | suffix-exact-match         | Items that end with `.js`              |
 * | `!.go$`     | inverse-suffix-exact-match | Items that do not end with `.go`       |
 *
 * A single pipe character acts as an OR operator. For example, the following
 * query matches entries that start with `core` and end with either`go`, `rb`,
 * or`py`.
 *
 * ```
 * ^core go$ | rb$ | py$
 * ```
 */
class ExtendedSearch {
  constructor(
    pattern,
    {
      isCaseSensitive = Config.isCaseSensitive,
      includeMatches = Config.includeMatches,
      minMatchCharLength = Config.minMatchCharLength,
      ignoreLocation = Config.ignoreLocation,
      findAllMatches = Config.findAllMatches,
      location = Config.location,
      threshold = Config.threshold,
      distance = Config.distance
    } = {}
  ) {
    this.query = null;
    this.options = {
      isCaseSensitive,
      includeMatches,
      minMatchCharLength,
      findAllMatches,
      ignoreLocation,
      location,
      threshold,
      distance
    };

    this.pattern = isCaseSensitive ? pattern : pattern.toLowerCase();
    this.query = parseQuery(this.pattern, this.options);
  }

  static condition(_, options) {
    return options.useExtendedSearch
  }

  searchIn(text) {
    const query = this.query;

    if (!query) {
      return {
        isMatch: false,
        score: 1
      }
    }

    const { includeMatches, isCaseSensitive } = this.options;

    text = isCaseSensitive ? text : text.toLowerCase();

    let numMatches = 0;
    let allIndices = [];
    let totalScore = 0;

    // ORs
    for (let i = 0, qLen = query.length; i < qLen; i += 1) {
      const searchers = query[i];

      // Reset indices
      allIndices.length = 0;
      numMatches = 0;

      // ANDs
      for (let j = 0, pLen = searchers.length; j < pLen; j += 1) {
        const searcher = searchers[j];
        const { isMatch, indices, score } = searcher.search(text);

        if (isMatch) {
          numMatches += 1;
          totalScore += score;
          if (includeMatches) {
            const type = searcher.constructor.type;
            if (MultiMatchSet.has(type)) {
              allIndices = [...allIndices, ...indices];
            } else {
              allIndices.push(indices);
            }
          }
        } else {
          totalScore = 0;
          numMatches = 0;
          allIndices.length = 0;
          break
        }
      }

      // OR condition, so if TRUE, return
      if (numMatches) {
        let result = {
          isMatch: true,
          score: totalScore / numMatches
        };

        if (includeMatches) {
          result.indices = allIndices;
        }

        return result
      }
    }

    // Nothing was matched
    return {
      isMatch: false,
      score: 1
    }
  }
}

const registeredSearchers = [];

function register(...args) {
  registeredSearchers.push(...args);
}

function createSearcher(pattern, options) {
  for (let i = 0, len = registeredSearchers.length; i < len; i += 1) {
    let searcherClass = registeredSearchers[i];
    if (searcherClass.condition(pattern, options)) {
      return new searcherClass(pattern, options)
    }
  }

  return new BitapSearch(pattern, options)
}

const LogicalOperator = {
  AND: '$and',
  OR: '$or'
};

const KeyType = {
  PATH: '$path',
  PATTERN: '$val'
};

const isExpression = (query) =>
  !!(query[LogicalOperator.AND] || query[LogicalOperator.OR]);

const isPath = (query) => !!query[KeyType.PATH];

const isLeaf = (query) =>
  !isArray(query) && isObject(query) && !isExpression(query);

const convertToExplicit = (query) => ({
  [LogicalOperator.AND]: Object.keys(query).map((key) => ({
    [key]: query[key]
  }))
});

// When `auto` is `true`, the parse function will infer and initialize and add
// the appropriate `Searcher` instance
function parse$1(query, options, { auto = true } = {}) {
  const next = (query) => {
    let keys = Object.keys(query);

    const isQueryPath = isPath(query);

    if (!isQueryPath && keys.length > 1 && !isExpression(query)) {
      return next(convertToExplicit(query))
    }

    if (isLeaf(query)) {
      const key = isQueryPath ? query[KeyType.PATH] : keys[0];

      const pattern = isQueryPath ? query[KeyType.PATTERN] : query[key];

      if (!isString(pattern)) {
        throw new Error(LOGICAL_SEARCH_INVALID_QUERY_FOR_KEY(key))
      }

      const obj = {
        keyId: createKeyId(key),
        pattern
      };

      if (auto) {
        obj.searcher = createSearcher(pattern, options);
      }

      return obj
    }

    let node = {
      children: [],
      operator: keys[0]
    };

    keys.forEach((key) => {
      const value = query[key];

      if (isArray(value)) {
        value.forEach((item) => {
          node.children.push(next(item));
        });
      }
    });

    return node
  };

  if (!isExpression(query)) {
    query = convertToExplicit(query);
  }

  return next(query)
}

// Practical scoring function
function computeScore(
  results,
  { ignoreFieldNorm = Config.ignoreFieldNorm }
) {
  results.forEach((result) => {
    let totalScore = 1;

    result.matches.forEach(({ key, norm, score }) => {
      const weight = key ? key.weight : null;

      totalScore *= Math.pow(
        score === 0 && weight ? Number.EPSILON : score,
        (weight || 1) * (ignoreFieldNorm ? 1 : norm)
      );
    });

    result.score = totalScore;
  });
}

function transformMatches(result, data) {
  const matches = result.matches;
  data.matches = [];

  if (!isDefined(matches)) {
    return
  }

  matches.forEach((match) => {
    if (!isDefined(match.indices) || !match.indices.length) {
      return
    }

    const { indices, value } = match;

    let obj = {
      indices,
      value
    };

    if (match.key) {
      obj.key = match.key.src;
    }

    if (match.idx > -1) {
      obj.refIndex = match.idx;
    }

    data.matches.push(obj);
  });
}

function transformScore(result, data) {
  data.score = result.score;
}

function format(
  results,
  docs,
  {
    includeMatches = Config.includeMatches,
    includeScore = Config.includeScore
  } = {}
) {
  const transformers = [];

  if (includeMatches) transformers.push(transformMatches);
  if (includeScore) transformers.push(transformScore);

  return results.map((result) => {
    const { idx } = result;

    const data = {
      item: docs[idx],
      refIndex: idx
    };

    if (transformers.length) {
      transformers.forEach((transformer) => {
        transformer(result, data);
      });
    }

    return data
  })
}

class Fuse {
  constructor(docs, options = {}, index) {
    this.options = { ...Config, ...options };

    if (
      this.options.useExtendedSearch &&
      !true
    ) {
      throw new Error(EXTENDED_SEARCH_UNAVAILABLE)
    }

    this._keyStore = new KeyStore(this.options.keys);

    this.setCollection(docs, index);
  }

  setCollection(docs, index) {
    this._docs = docs;

    if (index && !(index instanceof FuseIndex)) {
      throw new Error(INCORRECT_INDEX_TYPE)
    }

    this._myIndex =
      index ||
      createIndex(this.options.keys, this._docs, {
        getFn: this.options.getFn,
        fieldNormWeight: this.options.fieldNormWeight
      });
  }

  add(doc) {
    if (!isDefined(doc)) {
      return
    }

    this._docs.push(doc);
    this._myIndex.add(doc);
  }

  remove(predicate = (/* doc, idx */) => false) {
    const results = [];

    for (let i = 0, len = this._docs.length; i < len; i += 1) {
      const doc = this._docs[i];
      if (predicate(doc, i)) {
        this.removeAt(i);
        i -= 1;
        len -= 1;

        results.push(doc);
      }
    }

    return results
  }

  removeAt(idx) {
    this._docs.splice(idx, 1);
    this._myIndex.removeAt(idx);
  }

  getIndex() {
    return this._myIndex
  }

  search(query, { limit = -1 } = {}) {
    const {
      includeMatches,
      includeScore,
      shouldSort,
      sortFn,
      ignoreFieldNorm
    } = this.options;

    let results = isString(query)
      ? isString(this._docs[0])
        ? this._searchStringList(query)
        : this._searchObjectList(query)
      : this._searchLogical(query);

    computeScore(results, { ignoreFieldNorm });

    if (shouldSort) {
      results.sort(sortFn);
    }

    if (isNumber(limit) && limit > -1) {
      results = results.slice(0, limit);
    }

    return format(results, this._docs, {
      includeMatches,
      includeScore
    })
  }

  _searchStringList(query) {
    const searcher = createSearcher(query, this.options);
    const { records } = this._myIndex;
    const results = [];

    // Iterate over every string in the index
    records.forEach(({ v: text, i: idx, n: norm }) => {
      if (!isDefined(text)) {
        return
      }

      const { isMatch, score, indices } = searcher.searchIn(text);

      if (isMatch) {
        results.push({
          item: text,
          idx,
          matches: [{ score, value: text, norm, indices }]
        });
      }
    });

    return results
  }

  _searchLogical(query) {

    const expression = parse$1(query, this.options);

    const evaluate = (node, item, idx) => {
      if (!node.children) {
        const { keyId, searcher } = node;

        const matches = this._findMatches({
          key: this._keyStore.get(keyId),
          value: this._myIndex.getValueForItemAtKeyId(item, keyId),
          searcher
        });

        if (matches && matches.length) {
          return [
            {
              idx,
              item,
              matches
            }
          ]
        }

        return []
      }

      const res = [];
      for (let i = 0, len = node.children.length; i < len; i += 1) {
        const child = node.children[i];
        const result = evaluate(child, item, idx);
        if (result.length) {
          res.push(...result);
        } else if (node.operator === LogicalOperator.AND) {
          return []
        }
      }
      return res
    };

    const records = this._myIndex.records;
    const resultMap = {};
    const results = [];

    records.forEach(({ $: item, i: idx }) => {
      if (isDefined(item)) {
        let expResults = evaluate(expression, item, idx);

        if (expResults.length) {
          // Dedupe when adding
          if (!resultMap[idx]) {
            resultMap[idx] = { idx, item, matches: [] };
            results.push(resultMap[idx]);
          }
          expResults.forEach(({ matches }) => {
            resultMap[idx].matches.push(...matches);
          });
        }
      }
    });

    return results
  }

  _searchObjectList(query) {
    const searcher = createSearcher(query, this.options);
    const { keys, records } = this._myIndex;
    const results = [];

    // List is Array<Object>
    records.forEach(({ $: item, i: idx }) => {
      if (!isDefined(item)) {
        return
      }

      let matches = [];

      // Iterate over every key (i.e, path), and fetch the value at that key
      keys.forEach((key, keyIndex) => {
        matches.push(
          ...this._findMatches({
            key,
            value: item[keyIndex],
            searcher
          })
        );
      });

      if (matches.length) {
        results.push({
          idx,
          item,
          matches
        });
      }
    });

    return results
  }
  _findMatches({ key, value, searcher }) {
    if (!isDefined(value)) {
      return []
    }

    let matches = [];

    if (isArray(value)) {
      value.forEach(({ v: text, i: idx, n: norm }) => {
        if (!isDefined(text)) {
          return
        }

        const { isMatch, score, indices } = searcher.searchIn(text);

        if (isMatch) {
          matches.push({
            score,
            key,
            value: text,
            idx,
            norm,
            indices
          });
        }
      });
    } else {
      const { v: text, n: norm } = value;

      const { isMatch, score, indices } = searcher.searchIn(text);

      if (isMatch) {
        matches.push({ score, key, value: text, norm, indices });
      }
    }

    return matches
  }
}

Fuse.version = '7.0.0';
Fuse.createIndex = createIndex;
Fuse.parseIndex = parseIndex;
Fuse.config = Config;

{
  Fuse.parseQuery = parse$1;
}

{
  register(ExtendedSearch);
}

var commonjsGlobal = typeof globalThis !== 'undefined' ? globalThis : typeof window !== 'undefined' ? window : typeof global !== 'undefined' ? global : typeof self !== 'undefined' ? self : {};

function getDefaultExportFromCjs (x) {
	return x && x.__esModule && Object.prototype.hasOwnProperty.call(x, 'default') ? x['default'] : x;
}

var src = {
	compareTwoStrings:compareTwoStrings,
	findBestMatch:findBestMatch
};

function compareTwoStrings(first, second) {
	first = first.replace(/\s+/g, '');
	second = second.replace(/\s+/g, '');

	if (first === second) return 1; // identical or empty
	if (first.length < 2 || second.length < 2) return 0; // if either is a 0-letter or 1-letter string

	let firstBigrams = new Map();
	for (let i = 0; i < first.length - 1; i++) {
		const bigram = first.substring(i, i + 2);
		const count = firstBigrams.has(bigram)
			? firstBigrams.get(bigram) + 1
			: 1;

		firstBigrams.set(bigram, count);
	}
	let intersectionSize = 0;
	for (let i = 0; i < second.length - 1; i++) {
		const bigram = second.substring(i, i + 2);
		const count = firstBigrams.has(bigram)
			? firstBigrams.get(bigram)
			: 0;

		if (count > 0) {
			firstBigrams.set(bigram, count - 1);
			intersectionSize++;
		}
	}

	return (2.0 * intersectionSize) / (first.length + second.length - 2);
}

function findBestMatch(mainString, targetStrings) {
	if (!areArgsValid(mainString, targetStrings)) throw new Error('Bad arguments: First argument should be a string, second should be an array of strings');
	
	const ratings = [];
	let bestMatchIndex = 0;

	for (let i = 0; i < targetStrings.length; i++) {
		const currentTargetString = targetStrings[i];
		const currentRating = compareTwoStrings(mainString, currentTargetString);
		ratings.push({target: currentTargetString, rating: currentRating});
		if (currentRating > ratings[bestMatchIndex].rating) {
			bestMatchIndex = i;
		}
	}
	
	
	const bestMatch = ratings[bestMatchIndex];
	
	return { ratings: ratings, bestMatch: bestMatch, bestMatchIndex: bestMatchIndex };
}

function areArgsValid(mainString, targetStrings) {
	if (typeof mainString !== 'string') return false;
	if (!Array.isArray(targetStrings)) return false;
	if (!targetStrings.length) return false;
	if (targetStrings.find( function (s) { return typeof s !== 'string'})) return false;
	return true;
}

var stringSimilarity = /*@__PURE__*/getDefaultExportFromCjs(src);

const prepareFuseString = (input) => {
    input = input + " ";
    input = input.toLowerCase();
    input = input.replace(/\s\|\s/g, " ");
    input = input.replace(/(^|\s)[!']\S+/g, (match) => match.replace(/[!']/g, ""));
    input = input.replace(/(^|\s)-\S+/g, (match) => match.replace("-", "!"));
    input = moveNegativeMatchesToFront(input);
    return input;
};
const addQuoteIfNeeded = (input) => {
    return input.replace(/(^|[^-])"([^"]+)"/g, (match, p1, p2) => {
        return `${p1} '"${p2}"`;
    });
};
const moveNegativeMatchesToFront = (input) => {
    const negativePhrasePattern = /\s-!"[^"]*"\s*/g;
    const negativePhrases = input.match(negativePhrasePattern) || [];
    let modifiedInput = input.replace(negativePhrasePattern, " ");
    const negativeWordsPattern = /\s![^ ]*\s*/g;
    const negativeWords = modifiedInput.match(negativeWordsPattern) || [];
    modifiedInput = modifiedInput.replace(negativeWordsPattern, " ");
    modifiedInput = addQuoteIfNeeded(modifiedInput);
    const movedToFront = [...negativePhrases, ...negativeWords].join(" ") + modifiedInput;
    return movedToFront;
};
const extractWords = (str) => {
    const phrasePattern = /"[^"]*"|\S+/g;
    const words = [];
    let phrase = phrasePattern.exec(str);
    while (phrase) {
        words.push(phrase[0]);
        phrase = phrasePattern.exec(str);
    }
    return words;
};

class FuseSearcher {
    constructor(_indexDataProvider) {
        this._indexDataProvider = _indexDataProvider;
        this._fuseSearchConfig = {
            nameSimilarityOffset: 0.8,
            paragraphSimilarityOffset: 0.8,
            paragraphOffset: 30,
        };
        this._excludedWordsCount = 0;
        this._fuses = {};
    }
    resetAllCatalogs() {
        return __awaiter(this, void 0, void 0, function* () {
            yield this._indexDataProvider.deleteCatalogs();
            this._fuses = {};
        });
    }
    searchAll(query, ids) {
        return __awaiter(this, void 0, void 0, function* () {
            let result = [];
            for (const catalogName in ids)
                result = result.concat(yield this.search(query, catalogName, ids[catalogName]));
            return result.sort((a, b) => {
                if (a.count === b.count)
                    return a.score - b.score;
                return b.count - a.count;
            });
        });
    }
    search(query, catalogName, articleIds) {
        return __awaiter(this, void 0, void 0, function* () {
            var _a, _b;
            const quotesCount = (query.match(/"/g) || []).length;
            if (quotesCount % 2 !== 0) {
                const lastQuoteIndex = query.lastIndexOf('"');
                if (lastQuoteIndex !== -1) {
                    query = query.substring(0, lastQuoteIndex) + query.substring(lastQuoteIndex + 1);
                }
            }
            const doesQueryHaveOnlyNegativeWords = ((_a = query.match(/-"[^"]*"\s*/g)) === null || _a === void 0 ? void 0 : _a.join("")) === query || ((_b = query.match(/-[^ ]*\s*/g)) === null || _b === void 0 ? void 0 : _b.join("")) === query;
            if (!(articleIds === null || articleIds === void 0 ? void 0 : articleIds.length) || doesQueryHaveOnlyNegativeWords)
                return [];
            if (!this._fuses[catalogName])
                yield this._initFuses(catalogName);
            this._excludedWordsCount = 0;
            query = prepareFuseString(query);
            const searchedArray = this._fuses[catalogName].search(query);
            query = this._removeNegatedPhrases(query);
            query = this._removeExtraCharacters(query);
            const result = searchedArray
                .filter((searched) => articleIds.some((id) => id === searched.item.path))
                .map((searched) => this._processSearchResult(searched, query))
                .filter((searched) => { var _a, _b; return ((_a = searched.paragraph) === null || _a === void 0 ? void 0 : _a.length) !== 0 || ((_b = searched.name.targets) === null || _b === void 0 ? void 0 : _b.length) !== 0; });
            return result.sort((a, b) => {
                if (a.count === b.count)
                    return a.score - b.score;
                return b.count - a.count;
            });
        });
    }
    _initFuses(catalogName) {
        return __awaiter(this, void 0, void 0, function* () {
            this._fuses[catalogName] = new Fuse(yield this._indexDataProvider.getCatalogValue(catalogName), {
                keys: [{ name: "path" }, { name: "logicPath" }, { name: "title" }, { name: "content" }, { name: "tags" }],
                useExtendedSearch: true,
                includeScore: true,
                includeMatches: true,
                ignoreLocation: true,
                ignoreFieldNorm: true,
                threshold: 0.2,
                shouldSort: true,
            });
        });
    }
    _processSearchResult(searched, query) {
        var _a;
        const searchItem = {
            name: {},
            paragraph: [],
            count: 0,
            score: searched.score,
            url: (_a = searched.item.pathname) !== null && _a !== void 0 ? _a : "",
        };
        searchItem.name = this._getName(searched, query);
        searchItem.paragraph = this._getMatchingParagraphs(searched, query);
        searchItem.count = searchItem.name.targets.length + searchItem.paragraph.length;
        return searchItem;
    }
    _getName(searched, query) {
        var _a, _b;
        const title = (_b = (_a = searched.item.title) === null || _a === void 0 ? void 0 : _a.toString()) !== null && _b !== void 0 ? _b : "";
        const endNameIndex = title.length;
        const name = { targets: [], end: "" };
        const queryArray = this._getStringArray(query);
        let startNameIndex = 0;
        let endTargetIndex = 0;
        searched.matches.forEach((match) => {
            if (match.key === "title") {
                const merged = [];
                const currentInterval = this._getCurrentInterval(match.indices, merged);
                merged.push(currentInterval);
                for (let i = 0; i < merged.length; i++) {
                    const beginTargetIndex = merged[i][0];
                    const endTargetIndexRaw = merged[i][1];
                    const tempEndTargetIndex = endTargetIndexRaw + 1;
                    const start = title.slice(startNameIndex, beginTargetIndex);
                    const target = title.slice(beginTargetIndex, tempEndTargetIndex);
                    for (let j = 0; j < queryArray.length; j++) {
                        const similarity = stringSimilarity.compareTwoStrings(queryArray[j], target.toLowerCase());
                        if (similarity > this._fuseSearchConfig.nameSimilarityOffset) {
                            startNameIndex = tempEndTargetIndex;
                            endTargetIndex = tempEndTargetIndex;
                            name.targets.push({ start, target });
                            merged.splice(i, 1);
                            i--;
                            break;
                        }
                    }
                }
            }
        });
        name.end = title.slice(endTargetIndex, endNameIndex);
        return name;
    }
    _getCurrentInterval(indices, merged) {
        const tempIndices = [...indices];
        for (let i = 0; i < this._excludedWordsCount; i++)
            tempIndices.shift();
        tempIndices.sort((a, b) => a[0] - b[0]);
        let currentInterval = tempIndices[0];
        for (let i = 1; i < tempIndices.length; i++) {
            const nextInterval = tempIndices[i];
            if (currentInterval[1] >= nextInterval[0]) {
                currentInterval = [
                    Math.min(currentInterval[0], nextInterval[0]),
                    Math.max(currentInterval[1], nextInterval[1]),
                ];
            }
            else {
                merged.push(currentInterval);
                currentInterval = nextInterval;
            }
        }
        return currentInterval;
    }
    _getMatchingParagraphs(searched, query) {
        var _a;
        const content = (_a = searched.item.content) !== null && _a !== void 0 ? _a : "";
        const paragraphs = [];
        const queryArray = this._getStringArray(query);
        searched.matches.forEach((match) => {
            const merged = [];
            const currentInterval = this._getCurrentInterval(match.indices, merged);
            merged.push(currentInterval);
            if (match.key === "content") {
                merged.forEach((matchIndex) => {
                    for (let i = 0; i < queryArray.length; i++) {
                        const paragraph = this._getMatchingParagraph(content, matchIndex, queryArray[i]);
                        if (paragraph) {
                            paragraphs.push(paragraph);
                            break;
                        }
                    }
                });
            }
        });
        return paragraphs;
    }
    _getMatchingParagraph(text, matchIndex, query) {
        const [startIndex, endIndex] = matchIndex;
        const prevStart = Math.max(0, startIndex - this._fuseSearchConfig.paragraphOffset);
        const nextEnd = Math.min(text.length - 1, endIndex + this._fuseSearchConfig.paragraphOffset);
        const prev = prevStart > 0 ? "..." + text.slice(prevStart, startIndex) : text.slice(prevStart, startIndex);
        const target = text.slice(startIndex, endIndex + 1);
        const next = endIndex < nextEnd ? text.slice(endIndex + 1, nextEnd + 1) + "..." : text.slice(endIndex + 1, nextEnd + 1);
        const similarity = stringSimilarity.compareTwoStrings(query, target.toLowerCase());
        return similarity > this._fuseSearchConfig.paragraphSimilarityOffset ? { prev, target, next } : null;
    }
    _getStringArray(str) {
        const words = extractWords(str);
        const uniqueWords = [...new Set(words.map((word) => word.replace(/["']/g, "")))];
        return uniqueWords;
    }
    _removeNegatedPhrases(str) {
        var _a;
        const negativePhrasePattern = /!"[^"]*"\s*/g;
        const negativeWordsPattern = /![^ ]*\s*/g;
        const matchesPattern = str.match(negativeWordsPattern);
        this._excludedWordsCount += (_a = matchesPattern === null || matchesPattern === void 0 ? void 0 : matchesPattern.length) !== null && _a !== void 0 ? _a : 0;
        return str.replace(negativePhrasePattern, "").replace(negativeWordsPattern, "");
    }
    _removeExtraCharacters(str) {
        return str.replace(/ (['=^])/g, " ").replace(/([$]) /g, " ");
    }
}

var hp2Builder$2 = {};

var lib$5 = {};

var lib$4 = {};

(function (exports) {
	Object.defineProperty(exports, "__esModule", { value: true });
	exports.Doctype = exports.CDATA = exports.Tag = exports.Style = exports.Script = exports.Comment = exports.Directive = exports.Text = exports.Root = exports.isTag = exports.ElementType = void 0;
	/** Types of elements found in htmlparser2's DOM */
	var ElementType;
	(function (ElementType) {
	    /** Type for the root element of a document */
	    ElementType["Root"] = "root";
	    /** Type for Text */
	    ElementType["Text"] = "text";
	    /** Type for <? ... ?> */
	    ElementType["Directive"] = "directive";
	    /** Type for <!-- ... --> */
	    ElementType["Comment"] = "comment";
	    /** Type for <script> tags */
	    ElementType["Script"] = "script";
	    /** Type for <style> tags */
	    ElementType["Style"] = "style";
	    /** Type for Any tag */
	    ElementType["Tag"] = "tag";
	    /** Type for <![CDATA[ ... ]]> */
	    ElementType["CDATA"] = "cdata";
	    /** Type for <!doctype ...> */
	    ElementType["Doctype"] = "doctype";
	})(ElementType = exports.ElementType || (exports.ElementType = {}));
	/**
	 * Tests whether an element is a tag or not.
	 *
	 * @param elem Element to test
	 */
	function isTag(elem) {
	    return (elem.type === ElementType.Tag ||
	        elem.type === ElementType.Script ||
	        elem.type === ElementType.Style);
	}
	exports.isTag = isTag;
	// Exports for backwards compatibility
	/** Type for the root element of a document */
	exports.Root = ElementType.Root;
	/** Type for Text */
	exports.Text = ElementType.Text;
	/** Type for <? ... ?> */
	exports.Directive = ElementType.Directive;
	/** Type for <!-- ... --> */
	exports.Comment = ElementType.Comment;
	/** Type for <script> tags */
	exports.Script = ElementType.Script;
	/** Type for <style> tags */
	exports.Style = ElementType.Style;
	/** Type for Any tag */
	exports.Tag = ElementType.Tag;
	/** Type for <![CDATA[ ... ]]> */
	exports.CDATA = ElementType.CDATA;
	/** Type for <!doctype ...> */
	exports.Doctype = ElementType.Doctype; 
} (lib$4));

var node = {};

var __extends$1 = (commonjsGlobal && commonjsGlobal.__extends) || (function () {
    var extendStatics = function (d, b) {
        extendStatics = Object.setPrototypeOf ||
            ({ __proto__: [] } instanceof Array && function (d, b) { d.__proto__ = b; }) ||
            function (d, b) { for (var p in b) if (Object.prototype.hasOwnProperty.call(b, p)) d[p] = b[p]; };
        return extendStatics(d, b);
    };
    return function (d, b) {
        if (typeof b !== "function" && b !== null)
            throw new TypeError("Class extends value " + String(b) + " is not a constructor or null");
        extendStatics(d, b);
        function __() { this.constructor = d; }
        d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
    };
})();
var __assign$1 = (commonjsGlobal && commonjsGlobal.__assign) || function () {
    __assign$1 = Object.assign || function(t) {
        for (var s, i = 1, n = arguments.length; i < n; i++) {
            s = arguments[i];
            for (var p in s) if (Object.prototype.hasOwnProperty.call(s, p))
                t[p] = s[p];
        }
        return t;
    };
    return __assign$1.apply(this, arguments);
};
Object.defineProperty(node, "__esModule", { value: true });
node.cloneNode = node.hasChildren = node.isDocument = node.isDirective = node.isComment = node.isText = node.isCDATA = node.isTag = node.Element = node.Document = node.NodeWithChildren = node.ProcessingInstruction = node.Comment = node.Text = node.DataNode = node.Node = void 0;
var domelementtype_1$1 = lib$4;
var nodeTypes = new Map([
    [domelementtype_1$1.ElementType.Tag, 1],
    [domelementtype_1$1.ElementType.Script, 1],
    [domelementtype_1$1.ElementType.Style, 1],
    [domelementtype_1$1.ElementType.Directive, 1],
    [domelementtype_1$1.ElementType.Text, 3],
    [domelementtype_1$1.ElementType.CDATA, 4],
    [domelementtype_1$1.ElementType.Comment, 8],
    [domelementtype_1$1.ElementType.Root, 9],
]);
/**
 * This object will be used as the prototype for Nodes when creating a
 * DOM-Level-1-compliant structure.
 */
var Node = /** @class */ (function () {
    /**
     *
     * @param type The type of the node.
     */
    function Node(type) {
        this.type = type;
        /** Parent of the node */
        this.parent = null;
        /** Previous sibling */
        this.prev = null;
        /** Next sibling */
        this.next = null;
        /** The start index of the node. Requires `withStartIndices` on the handler to be `true. */
        this.startIndex = null;
        /** The end index of the node. Requires `withEndIndices` on the handler to be `true. */
        this.endIndex = null;
    }
    Object.defineProperty(Node.prototype, "nodeType", {
        // Read-only aliases
        /**
         * [DOM spec](https://dom.spec.whatwg.org/#dom-node-nodetype)-compatible
         * node {@link type}.
         */
        get: function () {
            var _a;
            return (_a = nodeTypes.get(this.type)) !== null && _a !== void 0 ? _a : 1;
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(Node.prototype, "parentNode", {
        // Read-write aliases for properties
        /**
         * Same as {@link parent}.
         * [DOM spec](https://dom.spec.whatwg.org)-compatible alias.
         */
        get: function () {
            return this.parent;
        },
        set: function (parent) {
            this.parent = parent;
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(Node.prototype, "previousSibling", {
        /**
         * Same as {@link prev}.
         * [DOM spec](https://dom.spec.whatwg.org)-compatible alias.
         */
        get: function () {
            return this.prev;
        },
        set: function (prev) {
            this.prev = prev;
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(Node.prototype, "nextSibling", {
        /**
         * Same as {@link next}.
         * [DOM spec](https://dom.spec.whatwg.org)-compatible alias.
         */
        get: function () {
            return this.next;
        },
        set: function (next) {
            this.next = next;
        },
        enumerable: false,
        configurable: true
    });
    /**
     * Clone this node, and optionally its children.
     *
     * @param recursive Clone child nodes as well.
     * @returns A clone of the node.
     */
    Node.prototype.cloneNode = function (recursive) {
        if (recursive === void 0) { recursive = false; }
        return cloneNode(this, recursive);
    };
    return Node;
}());
node.Node = Node;
/**
 * A node that contains some data.
 */
var DataNode = /** @class */ (function (_super) {
    __extends$1(DataNode, _super);
    /**
     * @param type The type of the node
     * @param data The content of the data node
     */
    function DataNode(type, data) {
        var _this = _super.call(this, type) || this;
        _this.data = data;
        return _this;
    }
    Object.defineProperty(DataNode.prototype, "nodeValue", {
        /**
         * Same as {@link data}.
         * [DOM spec](https://dom.spec.whatwg.org)-compatible alias.
         */
        get: function () {
            return this.data;
        },
        set: function (data) {
            this.data = data;
        },
        enumerable: false,
        configurable: true
    });
    return DataNode;
}(Node));
node.DataNode = DataNode;
/**
 * Text within the document.
 */
var Text = /** @class */ (function (_super) {
    __extends$1(Text, _super);
    function Text(data) {
        return _super.call(this, domelementtype_1$1.ElementType.Text, data) || this;
    }
    return Text;
}(DataNode));
node.Text = Text;
/**
 * Comments within the document.
 */
var Comment = /** @class */ (function (_super) {
    __extends$1(Comment, _super);
    function Comment(data) {
        return _super.call(this, domelementtype_1$1.ElementType.Comment, data) || this;
    }
    return Comment;
}(DataNode));
node.Comment = Comment;
/**
 * Processing instructions, including doc types.
 */
var ProcessingInstruction = /** @class */ (function (_super) {
    __extends$1(ProcessingInstruction, _super);
    function ProcessingInstruction(name, data) {
        var _this = _super.call(this, domelementtype_1$1.ElementType.Directive, data) || this;
        _this.name = name;
        return _this;
    }
    return ProcessingInstruction;
}(DataNode));
node.ProcessingInstruction = ProcessingInstruction;
/**
 * A `Node` that can have children.
 */
var NodeWithChildren = /** @class */ (function (_super) {
    __extends$1(NodeWithChildren, _super);
    /**
     * @param type Type of the node.
     * @param children Children of the node. Only certain node types can have children.
     */
    function NodeWithChildren(type, children) {
        var _this = _super.call(this, type) || this;
        _this.children = children;
        return _this;
    }
    Object.defineProperty(NodeWithChildren.prototype, "firstChild", {
        // Aliases
        /** First child of the node. */
        get: function () {
            var _a;
            return (_a = this.children[0]) !== null && _a !== void 0 ? _a : null;
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(NodeWithChildren.prototype, "lastChild", {
        /** Last child of the node. */
        get: function () {
            return this.children.length > 0
                ? this.children[this.children.length - 1]
                : null;
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(NodeWithChildren.prototype, "childNodes", {
        /**
         * Same as {@link children}.
         * [DOM spec](https://dom.spec.whatwg.org)-compatible alias.
         */
        get: function () {
            return this.children;
        },
        set: function (children) {
            this.children = children;
        },
        enumerable: false,
        configurable: true
    });
    return NodeWithChildren;
}(Node));
node.NodeWithChildren = NodeWithChildren;
/**
 * The root node of the document.
 */
var Document = /** @class */ (function (_super) {
    __extends$1(Document, _super);
    function Document(children) {
        return _super.call(this, domelementtype_1$1.ElementType.Root, children) || this;
    }
    return Document;
}(NodeWithChildren));
node.Document = Document;
/**
 * An element within the DOM.
 */
var Element$2 = /** @class */ (function (_super) {
    __extends$1(Element, _super);
    /**
     * @param name Name of the tag, eg. `div`, `span`.
     * @param attribs Object mapping attribute names to attribute values.
     * @param children Children of the node.
     */
    function Element(name, attribs, children, type) {
        if (children === void 0) { children = []; }
        if (type === void 0) { type = name === "script"
            ? domelementtype_1$1.ElementType.Script
            : name === "style"
                ? domelementtype_1$1.ElementType.Style
                : domelementtype_1$1.ElementType.Tag; }
        var _this = _super.call(this, type, children) || this;
        _this.name = name;
        _this.attribs = attribs;
        return _this;
    }
    Object.defineProperty(Element.prototype, "tagName", {
        // DOM Level 1 aliases
        /**
         * Same as {@link name}.
         * [DOM spec](https://dom.spec.whatwg.org)-compatible alias.
         */
        get: function () {
            return this.name;
        },
        set: function (name) {
            this.name = name;
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(Element.prototype, "attributes", {
        get: function () {
            var _this = this;
            return Object.keys(this.attribs).map(function (name) {
                var _a, _b;
                return ({
                    name: name,
                    value: _this.attribs[name],
                    namespace: (_a = _this["x-attribsNamespace"]) === null || _a === void 0 ? void 0 : _a[name],
                    prefix: (_b = _this["x-attribsPrefix"]) === null || _b === void 0 ? void 0 : _b[name],
                });
            });
        },
        enumerable: false,
        configurable: true
    });
    return Element;
}(NodeWithChildren));
node.Element = Element$2;
/**
 * @param node Node to check.
 * @returns `true` if the node is a `Element`, `false` otherwise.
 */
function isTag(node) {
    return (0, domelementtype_1$1.isTag)(node);
}
node.isTag = isTag;
/**
 * @param node Node to check.
 * @returns `true` if the node has the type `CDATA`, `false` otherwise.
 */
function isCDATA(node) {
    return node.type === domelementtype_1$1.ElementType.CDATA;
}
node.isCDATA = isCDATA;
/**
 * @param node Node to check.
 * @returns `true` if the node has the type `Text`, `false` otherwise.
 */
function isText(node) {
    return node.type === domelementtype_1$1.ElementType.Text;
}
node.isText = isText;
/**
 * @param node Node to check.
 * @returns `true` if the node has the type `Comment`, `false` otherwise.
 */
function isComment(node) {
    return node.type === domelementtype_1$1.ElementType.Comment;
}
node.isComment = isComment;
/**
 * @param node Node to check.
 * @returns `true` if the node has the type `ProcessingInstruction`, `false` otherwise.
 */
function isDirective(node) {
    return node.type === domelementtype_1$1.ElementType.Directive;
}
node.isDirective = isDirective;
/**
 * @param node Node to check.
 * @returns `true` if the node has the type `ProcessingInstruction`, `false` otherwise.
 */
function isDocument(node) {
    return node.type === domelementtype_1$1.ElementType.Root;
}
node.isDocument = isDocument;
/**
 * @param node Node to check.
 * @returns `true` if the node is a `NodeWithChildren` (has children), `false` otherwise.
 */
function hasChildren(node) {
    return Object.prototype.hasOwnProperty.call(node, "children");
}
node.hasChildren = hasChildren;
/**
 * Clone a node, and optionally its children.
 *
 * @param recursive Clone child nodes as well.
 * @returns A clone of the node.
 */
function cloneNode(node, recursive) {
    if (recursive === void 0) { recursive = false; }
    var result;
    if (isText(node)) {
        result = new Text(node.data);
    }
    else if (isComment(node)) {
        result = new Comment(node.data);
    }
    else if (isTag(node)) {
        var children = recursive ? cloneChildren(node.children) : [];
        var clone_1 = new Element$2(node.name, __assign$1({}, node.attribs), children);
        children.forEach(function (child) { return (child.parent = clone_1); });
        if (node.namespace != null) {
            clone_1.namespace = node.namespace;
        }
        if (node["x-attribsNamespace"]) {
            clone_1["x-attribsNamespace"] = __assign$1({}, node["x-attribsNamespace"]);
        }
        if (node["x-attribsPrefix"]) {
            clone_1["x-attribsPrefix"] = __assign$1({}, node["x-attribsPrefix"]);
        }
        result = clone_1;
    }
    else if (isCDATA(node)) {
        var children = recursive ? cloneChildren(node.children) : [];
        var clone_2 = new NodeWithChildren(domelementtype_1$1.ElementType.CDATA, children);
        children.forEach(function (child) { return (child.parent = clone_2); });
        result = clone_2;
    }
    else if (isDocument(node)) {
        var children = recursive ? cloneChildren(node.children) : [];
        var clone_3 = new Document(children);
        children.forEach(function (child) { return (child.parent = clone_3); });
        if (node["x-mode"]) {
            clone_3["x-mode"] = node["x-mode"];
        }
        result = clone_3;
    }
    else if (isDirective(node)) {
        var instruction = new ProcessingInstruction(node.name, node.data);
        if (node["x-name"] != null) {
            instruction["x-name"] = node["x-name"];
            instruction["x-publicId"] = node["x-publicId"];
            instruction["x-systemId"] = node["x-systemId"];
        }
        result = instruction;
    }
    else {
        throw new Error("Not implemented yet: ".concat(node.type));
    }
    result.startIndex = node.startIndex;
    result.endIndex = node.endIndex;
    if (node.sourceCodeLocation != null) {
        result.sourceCodeLocation = node.sourceCodeLocation;
    }
    return result;
}
node.cloneNode = cloneNode;
function cloneChildren(childs) {
    var children = childs.map(function (child) { return cloneNode(child, true); });
    for (var i = 1; i < children.length; i++) {
        children[i].prev = children[i - 1];
        children[i - 1].next = children[i];
    }
    return children;
}

(function (exports) {
	var __createBinding = (commonjsGlobal && commonjsGlobal.__createBinding) || (Object.create ? (function(o, m, k, k2) {
	    if (k2 === undefined) k2 = k;
	    var desc = Object.getOwnPropertyDescriptor(m, k);
	    if (!desc || ("get" in desc ? !m.__esModule : desc.writable || desc.configurable)) {
	      desc = { enumerable: true, get: function() { return m[k]; } };
	    }
	    Object.defineProperty(o, k2, desc);
	}) : (function(o, m, k, k2) {
	    if (k2 === undefined) k2 = k;
	    o[k2] = m[k];
	}));
	var __exportStar = (commonjsGlobal && commonjsGlobal.__exportStar) || function(m, exports) {
	    for (var p in m) if (p !== "default" && !Object.prototype.hasOwnProperty.call(exports, p)) __createBinding(exports, m, p);
	};
	Object.defineProperty(exports, "__esModule", { value: true });
	exports.DomHandler = void 0;
	var domelementtype_1 = lib$4;
	var node_1 = node;
	__exportStar(node, exports);
	var reWhitespace = /\s+/g;
	// Default options
	var defaultOpts = {
	    normalizeWhitespace: false,
	    withStartIndices: false,
	    withEndIndices: false,
	    xmlMode: false,
	};
	var DomHandler = /** @class */ (function () {
	    /**
	     * @param callback Called once parsing has completed.
	     * @param options Settings for the handler.
	     * @param elementCB Callback whenever a tag is closed.
	     */
	    function DomHandler(callback, options, elementCB) {
	        /** The elements of the DOM */
	        this.dom = [];
	        /** The root element for the DOM */
	        this.root = new node_1.Document(this.dom);
	        /** Indicated whether parsing has been completed. */
	        this.done = false;
	        /** Stack of open tags. */
	        this.tagStack = [this.root];
	        /** A data node that is still being written to. */
	        this.lastNode = null;
	        /** Reference to the parser instance. Used for location information. */
	        this.parser = null;
	        // Make it possible to skip arguments, for backwards-compatibility
	        if (typeof options === "function") {
	            elementCB = options;
	            options = defaultOpts;
	        }
	        if (typeof callback === "object") {
	            options = callback;
	            callback = undefined;
	        }
	        this.callback = callback !== null && callback !== void 0 ? callback : null;
	        this.options = options !== null && options !== void 0 ? options : defaultOpts;
	        this.elementCB = elementCB !== null && elementCB !== void 0 ? elementCB : null;
	    }
	    DomHandler.prototype.onparserinit = function (parser) {
	        this.parser = parser;
	    };
	    // Resets the handler back to starting state
	    DomHandler.prototype.onreset = function () {
	        this.dom = [];
	        this.root = new node_1.Document(this.dom);
	        this.done = false;
	        this.tagStack = [this.root];
	        this.lastNode = null;
	        this.parser = null;
	    };
	    // Signals the handler that parsing is done
	    DomHandler.prototype.onend = function () {
	        if (this.done)
	            return;
	        this.done = true;
	        this.parser = null;
	        this.handleCallback(null);
	    };
	    DomHandler.prototype.onerror = function (error) {
	        this.handleCallback(error);
	    };
	    DomHandler.prototype.onclosetag = function () {
	        this.lastNode = null;
	        var elem = this.tagStack.pop();
	        if (this.options.withEndIndices) {
	            elem.endIndex = this.parser.endIndex;
	        }
	        if (this.elementCB)
	            this.elementCB(elem);
	    };
	    DomHandler.prototype.onopentag = function (name, attribs) {
	        var type = this.options.xmlMode ? domelementtype_1.ElementType.Tag : undefined;
	        var element = new node_1.Element(name, attribs, undefined, type);
	        this.addNode(element);
	        this.tagStack.push(element);
	    };
	    DomHandler.prototype.ontext = function (data) {
	        var normalizeWhitespace = this.options.normalizeWhitespace;
	        var lastNode = this.lastNode;
	        if (lastNode && lastNode.type === domelementtype_1.ElementType.Text) {
	            if (normalizeWhitespace) {
	                lastNode.data = (lastNode.data + data).replace(reWhitespace, " ");
	            }
	            else {
	                lastNode.data += data;
	            }
	            if (this.options.withEndIndices) {
	                lastNode.endIndex = this.parser.endIndex;
	            }
	        }
	        else {
	            if (normalizeWhitespace) {
	                data = data.replace(reWhitespace, " ");
	            }
	            var node = new node_1.Text(data);
	            this.addNode(node);
	            this.lastNode = node;
	        }
	    };
	    DomHandler.prototype.oncomment = function (data) {
	        if (this.lastNode && this.lastNode.type === domelementtype_1.ElementType.Comment) {
	            this.lastNode.data += data;
	            return;
	        }
	        var node = new node_1.Comment(data);
	        this.addNode(node);
	        this.lastNode = node;
	    };
	    DomHandler.prototype.oncommentend = function () {
	        this.lastNode = null;
	    };
	    DomHandler.prototype.oncdatastart = function () {
	        var text = new node_1.Text("");
	        var node = new node_1.NodeWithChildren(domelementtype_1.ElementType.CDATA, [text]);
	        this.addNode(node);
	        text.parent = node;
	        this.lastNode = text;
	    };
	    DomHandler.prototype.oncdataend = function () {
	        this.lastNode = null;
	    };
	    DomHandler.prototype.onprocessinginstruction = function (name, data) {
	        var node = new node_1.ProcessingInstruction(name, data);
	        this.addNode(node);
	    };
	    DomHandler.prototype.handleCallback = function (error) {
	        if (typeof this.callback === "function") {
	            this.callback(error, this.dom);
	        }
	        else if (error) {
	            throw error;
	        }
	    };
	    DomHandler.prototype.addNode = function (node) {
	        var parent = this.tagStack[this.tagStack.length - 1];
	        var previousSibling = parent.children[parent.children.length - 1];
	        if (this.options.withStartIndices) {
	            node.startIndex = this.parser.startIndex;
	        }
	        if (this.options.withEndIndices) {
	            node.endIndex = this.parser.endIndex;
	        }
	        parent.children.push(node);
	        if (previousSibling) {
	            node.prev = previousSibling;
	            previousSibling.next = node;
	        }
	        node.parent = parent;
	        this.lastNode = null;
	    };
	    return DomHandler;
	}());
	exports.DomHandler = DomHandler;
	exports.default = DomHandler; 
} (lib$5));

var selderee$2 = {};

var parseley$1 = {};

var nearley$1 = {exports: {}};

(function (module) {
	(function(root, factory) {
	    if (module.exports) {
	        module.exports = factory();
	    } else {
	        root.nearley = factory();
	    }
	}(commonjsGlobal, function() {

	    function Rule(name, symbols, postprocess) {
	        this.id = ++Rule.highestId;
	        this.name = name;
	        this.symbols = symbols;        // a list of literal | regex class | nonterminal
	        this.postprocess = postprocess;
	        return this;
	    }
	    Rule.highestId = 0;

	    Rule.prototype.toString = function(withCursorAt) {
	        var symbolSequence = (typeof withCursorAt === "undefined")
	                             ? this.symbols.map(getSymbolShortDisplay).join(' ')
	                             : (   this.symbols.slice(0, withCursorAt).map(getSymbolShortDisplay).join(' ')
	                                 + " ● "
	                                 + this.symbols.slice(withCursorAt).map(getSymbolShortDisplay).join(' ')     );
	        return this.name + " → " + symbolSequence;
	    };


	    // a State is a rule at a position from a given starting point in the input stream (reference)
	    function State(rule, dot, reference, wantedBy) {
	        this.rule = rule;
	        this.dot = dot;
	        this.reference = reference;
	        this.data = [];
	        this.wantedBy = wantedBy;
	        this.isComplete = this.dot === rule.symbols.length;
	    }

	    State.prototype.toString = function() {
	        return "{" + this.rule.toString(this.dot) + "}, from: " + (this.reference || 0);
	    };

	    State.prototype.nextState = function(child) {
	        var state = new State(this.rule, this.dot + 1, this.reference, this.wantedBy);
	        state.left = this;
	        state.right = child;
	        if (state.isComplete) {
	            state.data = state.build();
	            // Having right set here will prevent the right state and its children
	            // form being garbage collected
	            state.right = undefined;
	        }
	        return state;
	    };

	    State.prototype.build = function() {
	        var children = [];
	        var node = this;
	        do {
	            children.push(node.right.data);
	            node = node.left;
	        } while (node.left);
	        children.reverse();
	        return children;
	    };

	    State.prototype.finish = function() {
	        if (this.rule.postprocess) {
	            this.data = this.rule.postprocess(this.data, this.reference, Parser.fail);
	        }
	    };


	    function Column(grammar, index) {
	        this.grammar = grammar;
	        this.index = index;
	        this.states = [];
	        this.wants = {}; // states indexed by the non-terminal they expect
	        this.scannable = []; // list of states that expect a token
	        this.completed = {}; // states that are nullable
	    }


	    Column.prototype.process = function(nextColumn) {
	        var states = this.states;
	        var wants = this.wants;
	        var completed = this.completed;

	        for (var w = 0; w < states.length; w++) { // nb. we push() during iteration
	            var state = states[w];

	            if (state.isComplete) {
	                state.finish();
	                if (state.data !== Parser.fail) {
	                    // complete
	                    var wantedBy = state.wantedBy;
	                    for (var i = wantedBy.length; i--; ) { // this line is hot
	                        var left = wantedBy[i];
	                        this.complete(left, state);
	                    }

	                    // special-case nullables
	                    if (state.reference === this.index) {
	                        // make sure future predictors of this rule get completed.
	                        var exp = state.rule.name;
	                        (this.completed[exp] = this.completed[exp] || []).push(state);
	                    }
	                }

	            } else {
	                // queue scannable states
	                var exp = state.rule.symbols[state.dot];
	                if (typeof exp !== 'string') {
	                    this.scannable.push(state);
	                    continue;
	                }

	                // predict
	                if (wants[exp]) {
	                    wants[exp].push(state);

	                    if (completed.hasOwnProperty(exp)) {
	                        var nulls = completed[exp];
	                        for (var i = 0; i < nulls.length; i++) {
	                            var right = nulls[i];
	                            this.complete(state, right);
	                        }
	                    }
	                } else {
	                    wants[exp] = [state];
	                    this.predict(exp);
	                }
	            }
	        }
	    };

	    Column.prototype.predict = function(exp) {
	        var rules = this.grammar.byName[exp] || [];

	        for (var i = 0; i < rules.length; i++) {
	            var r = rules[i];
	            var wantedBy = this.wants[exp];
	            var s = new State(r, 0, this.index, wantedBy);
	            this.states.push(s);
	        }
	    };

	    Column.prototype.complete = function(left, right) {
	        var copy = left.nextState(right);
	        this.states.push(copy);
	    };


	    function Grammar(rules, start) {
	        this.rules = rules;
	        this.start = start || this.rules[0].name;
	        var byName = this.byName = {};
	        this.rules.forEach(function(rule) {
	            if (!byName.hasOwnProperty(rule.name)) {
	                byName[rule.name] = [];
	            }
	            byName[rule.name].push(rule);
	        });
	    }

	    // So we can allow passing (rules, start) directly to Parser for backwards compatibility
	    Grammar.fromCompiled = function(rules, start) {
	        var lexer = rules.Lexer;
	        if (rules.ParserStart) {
	          start = rules.ParserStart;
	          rules = rules.ParserRules;
	        }
	        var rules = rules.map(function (r) { return (new Rule(r.name, r.symbols, r.postprocess)); });
	        var g = new Grammar(rules, start);
	        g.lexer = lexer; // nb. storing lexer on Grammar is iffy, but unavoidable
	        return g;
	    };


	    function StreamLexer() {
	      this.reset("");
	    }

	    StreamLexer.prototype.reset = function(data, state) {
	        this.buffer = data;
	        this.index = 0;
	        this.line = state ? state.line : 1;
	        this.lastLineBreak = state ? -state.col : 0;
	    };

	    StreamLexer.prototype.next = function() {
	        if (this.index < this.buffer.length) {
	            var ch = this.buffer[this.index++];
	            if (ch === '\n') {
	              this.line += 1;
	              this.lastLineBreak = this.index;
	            }
	            return {value: ch};
	        }
	    };

	    StreamLexer.prototype.save = function() {
	      return {
	        line: this.line,
	        col: this.index - this.lastLineBreak,
	      }
	    };

	    StreamLexer.prototype.formatError = function(token, message) {
	        // nb. this gets called after consuming the offending token,
	        // so the culprit is index-1
	        var buffer = this.buffer;
	        if (typeof buffer === 'string') {
	            var lines = buffer
	                .split("\n")
	                .slice(
	                    Math.max(0, this.line - 5), 
	                    this.line
	                );

	            var nextLineBreak = buffer.indexOf('\n', this.index);
	            if (nextLineBreak === -1) nextLineBreak = buffer.length;
	            var col = this.index - this.lastLineBreak;
	            var lastLineDigits = String(this.line).length;
	            message += " at line " + this.line + " col " + col + ":\n\n";
	            message += lines
	                .map(function(line, i) {
	                    return pad(this.line - lines.length + i + 1, lastLineDigits) + " " + line;
	                }, this)
	                .join("\n");
	            message += "\n" + pad("", lastLineDigits + col) + "^\n";
	            return message;
	        } else {
	            return message + " at index " + (this.index - 1);
	        }

	        function pad(n, length) {
	            var s = String(n);
	            return Array(length - s.length + 1).join(" ") + s;
	        }
	    };

	    function Parser(rules, start, options) {
	        if (rules instanceof Grammar) {
	            var grammar = rules;
	            var options = start;
	        } else {
	            var grammar = Grammar.fromCompiled(rules, start);
	        }
	        this.grammar = grammar;

	        // Read options
	        this.options = {
	            keepHistory: false,
	            lexer: grammar.lexer || new StreamLexer,
	        };
	        for (var key in (options || {})) {
	            this.options[key] = options[key];
	        }

	        // Setup lexer
	        this.lexer = this.options.lexer;
	        this.lexerState = undefined;

	        // Setup a table
	        var column = new Column(grammar, 0);
	        this.table = [column];

	        // I could be expecting anything.
	        column.wants[grammar.start] = [];
	        column.predict(grammar.start);
	        // TODO what if start rule is nullable?
	        column.process();
	        this.current = 0; // token index
	    }

	    // create a reserved token for indicating a parse fail
	    Parser.fail = {};

	    Parser.prototype.feed = function(chunk) {
	        var lexer = this.lexer;
	        lexer.reset(chunk, this.lexerState);

	        var token;
	        while (true) {
	            try {
	                token = lexer.next();
	                if (!token) {
	                    break;
	                }
	            } catch (e) {
	                // Create the next column so that the error reporter
	                // can display the correctly predicted states.
	                var nextColumn = new Column(this.grammar, this.current + 1);
	                this.table.push(nextColumn);
	                var err = new Error(this.reportLexerError(e));
	                err.offset = this.current;
	                err.token = e.token;
	                throw err;
	            }
	            // We add new states to table[current+1]
	            var column = this.table[this.current];

	            // GC unused states
	            if (!this.options.keepHistory) {
	                delete this.table[this.current - 1];
	            }

	            var n = this.current + 1;
	            var nextColumn = new Column(this.grammar, n);
	            this.table.push(nextColumn);

	            // Advance all tokens that expect the symbol
	            var literal = token.text !== undefined ? token.text : token.value;
	            var value = lexer.constructor === StreamLexer ? token.value : token;
	            var scannable = column.scannable;
	            for (var w = scannable.length; w--; ) {
	                var state = scannable[w];
	                var expect = state.rule.symbols[state.dot];
	                // Try to consume the token
	                // either regex or literal
	                if (expect.test ? expect.test(value) :
	                    expect.type ? expect.type === token.type
	                                : expect.literal === literal) {
	                    // Add it
	                    var next = state.nextState({data: value, token: token, isToken: true, reference: n - 1});
	                    nextColumn.states.push(next);
	                }
	            }

	            // Next, for each of the rules, we either
	            // (a) complete it, and try to see if the reference row expected that
	            //     rule
	            // (b) predict the next nonterminal it expects by adding that
	            //     nonterminal's start state
	            // To prevent duplication, we also keep track of rules we have already
	            // added

	            nextColumn.process();

	            // If needed, throw an error:
	            if (nextColumn.states.length === 0) {
	                // No states at all! This is not good.
	                var err = new Error(this.reportError(token));
	                err.offset = this.current;
	                err.token = token;
	                throw err;
	            }

	            // maybe save lexer state
	            if (this.options.keepHistory) {
	              column.lexerState = lexer.save();
	            }

	            this.current++;
	        }
	        if (column) {
	          this.lexerState = lexer.save();
	        }

	        // Incrementally keep track of results
	        this.results = this.finish();

	        // Allow chaining, for whatever it's worth
	        return this;
	    };

	    Parser.prototype.reportLexerError = function(lexerError) {
	        var tokenDisplay, lexerMessage;
	        // Planning to add a token property to moo's thrown error
	        // even on erroring tokens to be used in error display below
	        var token = lexerError.token;
	        if (token) {
	            tokenDisplay = "input " + JSON.stringify(token.text[0]) + " (lexer error)";
	            lexerMessage = this.lexer.formatError(token, "Syntax error");
	        } else {
	            tokenDisplay = "input (lexer error)";
	            lexerMessage = lexerError.message;
	        }
	        return this.reportErrorCommon(lexerMessage, tokenDisplay);
	    };

	    Parser.prototype.reportError = function(token) {
	        var tokenDisplay = (token.type ? token.type + " token: " : "") + JSON.stringify(token.value !== undefined ? token.value : token);
	        var lexerMessage = this.lexer.formatError(token, "Syntax error");
	        return this.reportErrorCommon(lexerMessage, tokenDisplay);
	    };

	    Parser.prototype.reportErrorCommon = function(lexerMessage, tokenDisplay) {
	        var lines = [];
	        lines.push(lexerMessage);
	        var lastColumnIndex = this.table.length - 2;
	        var lastColumn = this.table[lastColumnIndex];
	        var expectantStates = lastColumn.states
	            .filter(function(state) {
	                var nextSymbol = state.rule.symbols[state.dot];
	                return nextSymbol && typeof nextSymbol !== "string";
	            });

	        if (expectantStates.length === 0) {
	            lines.push('Unexpected ' + tokenDisplay + '. I did not expect any more input. Here is the state of my parse table:\n');
	            this.displayStateStack(lastColumn.states, lines);
	        } else {
	            lines.push('Unexpected ' + tokenDisplay + '. Instead, I was expecting to see one of the following:\n');
	            // Display a "state stack" for each expectant state
	            // - which shows you how this state came to be, step by step.
	            // If there is more than one derivation, we only display the first one.
	            var stateStacks = expectantStates
	                .map(function(state) {
	                    return this.buildFirstStateStack(state, []) || [state];
	                }, this);
	            // Display each state that is expecting a terminal symbol next.
	            stateStacks.forEach(function(stateStack) {
	                var state = stateStack[0];
	                var nextSymbol = state.rule.symbols[state.dot];
	                var symbolDisplay = this.getSymbolDisplay(nextSymbol);
	                lines.push('A ' + symbolDisplay + ' based on:');
	                this.displayStateStack(stateStack, lines);
	            }, this);
	        }
	        lines.push("");
	        return lines.join("\n");
	    };
	    
	    Parser.prototype.displayStateStack = function(stateStack, lines) {
	        var lastDisplay;
	        var sameDisplayCount = 0;
	        for (var j = 0; j < stateStack.length; j++) {
	            var state = stateStack[j];
	            var display = state.rule.toString(state.dot);
	            if (display === lastDisplay) {
	                sameDisplayCount++;
	            } else {
	                if (sameDisplayCount > 0) {
	                    lines.push('    ^ ' + sameDisplayCount + ' more lines identical to this');
	                }
	                sameDisplayCount = 0;
	                lines.push('    ' + display);
	            }
	            lastDisplay = display;
	        }
	    };

	    Parser.prototype.getSymbolDisplay = function(symbol) {
	        return getSymbolLongDisplay(symbol);
	    };

	    /*
	    Builds a the first state stack. You can think of a state stack as the call stack
	    of the recursive-descent parser which the Nearley parse algorithm simulates.
	    A state stack is represented as an array of state objects. Within a
	    state stack, the first item of the array will be the starting
	    state, with each successive item in the array going further back into history.

	    This function needs to be given a starting state and an empty array representing
	    the visited states, and it returns an single state stack.

	    */
	    Parser.prototype.buildFirstStateStack = function(state, visited) {
	        if (visited.indexOf(state) !== -1) {
	            // Found cycle, return null
	            // to eliminate this path from the results, because
	            // we don't know how to display it meaningfully
	            return null;
	        }
	        if (state.wantedBy.length === 0) {
	            return [state];
	        }
	        var prevState = state.wantedBy[0];
	        var childVisited = [state].concat(visited);
	        var childResult = this.buildFirstStateStack(prevState, childVisited);
	        if (childResult === null) {
	            return null;
	        }
	        return [state].concat(childResult);
	    };

	    Parser.prototype.save = function() {
	        var column = this.table[this.current];
	        column.lexerState = this.lexerState;
	        return column;
	    };

	    Parser.prototype.restore = function(column) {
	        var index = column.index;
	        this.current = index;
	        this.table[index] = column;
	        this.table.splice(index + 1);
	        this.lexerState = column.lexerState;

	        // Incrementally keep track of results
	        this.results = this.finish();
	    };

	    // nb. deprecated: use save/restore instead!
	    Parser.prototype.rewind = function(index) {
	        if (!this.options.keepHistory) {
	            throw new Error('set option `keepHistory` to enable rewinding')
	        }
	        // nb. recall column (table) indicies fall between token indicies.
	        //        col 0   --   token 0   --   col 1
	        this.restore(this.table[index]);
	    };

	    Parser.prototype.finish = function() {
	        // Return the possible parsings
	        var considerations = [];
	        var start = this.grammar.start;
	        var column = this.table[this.table.length - 1];
	        column.states.forEach(function (t) {
	            if (t.rule.name === start
	                    && t.dot === t.rule.symbols.length
	                    && t.reference === 0
	                    && t.data !== Parser.fail) {
	                considerations.push(t);
	            }
	        });
	        return considerations.map(function(c) {return c.data; });
	    };

	    function getSymbolLongDisplay(symbol) {
	        var type = typeof symbol;
	        if (type === "string") {
	            return symbol;
	        } else if (type === "object") {
	            if (symbol.literal) {
	                return JSON.stringify(symbol.literal);
	            } else if (symbol instanceof RegExp) {
	                return 'character matching ' + symbol;
	            } else if (symbol.type) {
	                return symbol.type + ' token';
	            } else if (symbol.test) {
	                return 'token matching ' + String(symbol.test);
	            } else {
	                throw new Error('Unknown symbol type: ' + symbol);
	            }
	        }
	    }

	    function getSymbolShortDisplay(symbol) {
	        var type = typeof symbol;
	        if (type === "string") {
	            return symbol;
	        } else if (type === "object") {
	            if (symbol.literal) {
	                return JSON.stringify(symbol.literal);
	            } else if (symbol instanceof RegExp) {
	                return symbol.toString();
	            } else if (symbol.type) {
	                return '%' + symbol.type;
	            } else if (symbol.test) {
	                return '<' + String(symbol.test) + '>';
	            } else {
	                throw new Error('Unknown symbol type: ' + symbol);
	            }
	        }
	    }

	    return {
	        Parser: Parser,
	        Grammar: Grammar,
	        Rule: Rule,
	    };

	})); 
} (nearley$1));

var nearleyExports = nearley$1.exports;

var moo$1 = {exports: {}};

(function (module) {
	(function(root, factory) {
	  if (module.exports) {
	    module.exports = factory();
	  } else {
	    root.moo = factory();
	  }
	}(commonjsGlobal, function() {

	  var hasOwnProperty = Object.prototype.hasOwnProperty;
	  var toString = Object.prototype.toString;
	  var hasSticky = typeof new RegExp().sticky === 'boolean';

	  /***************************************************************************/

	  function isRegExp(o) { return o && toString.call(o) === '[object RegExp]' }
	  function isObject(o) { return o && typeof o === 'object' && !isRegExp(o) && !Array.isArray(o) }

	  function reEscape(s) {
	    return s.replace(/[-\/\\^$*+?.()|[\]{}]/g, '\\$&')
	  }
	  function reGroups(s) {
	    var re = new RegExp('|' + s);
	    return re.exec('').length - 1
	  }
	  function reCapture(s) {
	    return '(' + s + ')'
	  }
	  function reUnion(regexps) {
	    if (!regexps.length) return '(?!)'
	    var source =  regexps.map(function(s) {
	      return "(?:" + s + ")"
	    }).join('|');
	    return "(?:" + source + ")"
	  }

	  function regexpOrLiteral(obj) {
	    if (typeof obj === 'string') {
	      return '(?:' + reEscape(obj) + ')'

	    } else if (isRegExp(obj)) {
	      // TODO: consider /u support
	      if (obj.ignoreCase) throw new Error('RegExp /i flag not allowed')
	      if (obj.global) throw new Error('RegExp /g flag is implied')
	      if (obj.sticky) throw new Error('RegExp /y flag is implied')
	      if (obj.multiline) throw new Error('RegExp /m flag is implied')
	      return obj.source

	    } else {
	      throw new Error('Not a pattern: ' + obj)
	    }
	  }

	  function pad(s, length) {
	    if (s.length > length) {
	      return s
	    }
	    return Array(length - s.length + 1).join(" ") + s
	  }

	  function lastNLines(string, numLines) {
	    var position = string.length;
	    var lineBreaks = 0;
	    while (true) {
	      var idx = string.lastIndexOf("\n", position - 1);
	      if (idx === -1) {
	        break;
	      } else {
	        lineBreaks++;
	      }
	      position = idx;
	      if (lineBreaks === numLines) {
	        break;
	      }
	      if (position === 0) {
	        break;
	      }
	    }
	    var startPosition = 
	      lineBreaks < numLines ?
	      0 : 
	      position + 1;
	    return string.substring(startPosition).split("\n")
	  }

	  function objectToRules(object) {
	    var keys = Object.getOwnPropertyNames(object);
	    var result = [];
	    for (var i = 0; i < keys.length; i++) {
	      var key = keys[i];
	      var thing = object[key];
	      var rules = [].concat(thing);
	      if (key === 'include') {
	        for (var j = 0; j < rules.length; j++) {
	          result.push({include: rules[j]});
	        }
	        continue
	      }
	      var match = [];
	      rules.forEach(function(rule) {
	        if (isObject(rule)) {
	          if (match.length) result.push(ruleOptions(key, match));
	          result.push(ruleOptions(key, rule));
	          match = [];
	        } else {
	          match.push(rule);
	        }
	      });
	      if (match.length) result.push(ruleOptions(key, match));
	    }
	    return result
	  }

	  function arrayToRules(array) {
	    var result = [];
	    for (var i = 0; i < array.length; i++) {
	      var obj = array[i];
	      if (obj.include) {
	        var include = [].concat(obj.include);
	        for (var j = 0; j < include.length; j++) {
	          result.push({include: include[j]});
	        }
	        continue
	      }
	      if (!obj.type) {
	        throw new Error('Rule has no type: ' + JSON.stringify(obj))
	      }
	      result.push(ruleOptions(obj.type, obj));
	    }
	    return result
	  }

	  function ruleOptions(type, obj) {
	    if (!isObject(obj)) {
	      obj = { match: obj };
	    }
	    if (obj.include) {
	      throw new Error('Matching rules cannot also include states')
	    }

	    // nb. error and fallback imply lineBreaks
	    var options = {
	      defaultType: type,
	      lineBreaks: !!obj.error || !!obj.fallback,
	      pop: false,
	      next: null,
	      push: null,
	      error: false,
	      fallback: false,
	      value: null,
	      type: null,
	      shouldThrow: false,
	    };

	    // Avoid Object.assign(), so we support IE9+
	    for (var key in obj) {
	      if (hasOwnProperty.call(obj, key)) {
	        options[key] = obj[key];
	      }
	    }

	    // type transform cannot be a string
	    if (typeof options.type === 'string' && type !== options.type) {
	      throw new Error("Type transform cannot be a string (type '" + options.type + "' for token '" + type + "')")
	    }

	    // convert to array
	    var match = options.match;
	    options.match = Array.isArray(match) ? match : match ? [match] : [];
	    options.match.sort(function(a, b) {
	      return isRegExp(a) && isRegExp(b) ? 0
	           : isRegExp(b) ? -1 : isRegExp(a) ? +1 : b.length - a.length
	    });
	    return options
	  }

	  function toRules(spec) {
	    return Array.isArray(spec) ? arrayToRules(spec) : objectToRules(spec)
	  }

	  var defaultErrorRule = ruleOptions('error', {lineBreaks: true, shouldThrow: true});
	  function compileRules(rules, hasStates) {
	    var errorRule = null;
	    var fast = Object.create(null);
	    var fastAllowed = true;
	    var unicodeFlag = null;
	    var groups = [];
	    var parts = [];

	    // If there is a fallback rule, then disable fast matching
	    for (var i = 0; i < rules.length; i++) {
	      if (rules[i].fallback) {
	        fastAllowed = false;
	      }
	    }

	    for (var i = 0; i < rules.length; i++) {
	      var options = rules[i];

	      if (options.include) {
	        // all valid inclusions are removed by states() preprocessor
	        throw new Error('Inheritance is not allowed in stateless lexers')
	      }

	      if (options.error || options.fallback) {
	        // errorRule can only be set once
	        if (errorRule) {
	          if (!options.fallback === !errorRule.fallback) {
	            throw new Error("Multiple " + (options.fallback ? "fallback" : "error") + " rules not allowed (for token '" + options.defaultType + "')")
	          } else {
	            throw new Error("fallback and error are mutually exclusive (for token '" + options.defaultType + "')")
	          }
	        }
	        errorRule = options;
	      }

	      var match = options.match.slice();
	      if (fastAllowed) {
	        while (match.length && typeof match[0] === 'string' && match[0].length === 1) {
	          var word = match.shift();
	          fast[word.charCodeAt(0)] = options;
	        }
	      }

	      // Warn about inappropriate state-switching options
	      if (options.pop || options.push || options.next) {
	        if (!hasStates) {
	          throw new Error("State-switching options are not allowed in stateless lexers (for token '" + options.defaultType + "')")
	        }
	        if (options.fallback) {
	          throw new Error("State-switching options are not allowed on fallback tokens (for token '" + options.defaultType + "')")
	        }
	      }

	      // Only rules with a .match are included in the RegExp
	      if (match.length === 0) {
	        continue
	      }
	      fastAllowed = false;

	      groups.push(options);

	      // Check unicode flag is used everywhere or nowhere
	      for (var j = 0; j < match.length; j++) {
	        var obj = match[j];
	        if (!isRegExp(obj)) {
	          continue
	        }

	        if (unicodeFlag === null) {
	          unicodeFlag = obj.unicode;
	        } else if (unicodeFlag !== obj.unicode && options.fallback === false) {
	          throw new Error('If one rule is /u then all must be')
	        }
	      }

	      // convert to RegExp
	      var pat = reUnion(match.map(regexpOrLiteral));

	      // validate
	      var regexp = new RegExp(pat);
	      if (regexp.test("")) {
	        throw new Error("RegExp matches empty string: " + regexp)
	      }
	      var groupCount = reGroups(pat);
	      if (groupCount > 0) {
	        throw new Error("RegExp has capture groups: " + regexp + "\nUse (?: … ) instead")
	      }

	      // try and detect rules matching newlines
	      if (!options.lineBreaks && regexp.test('\n')) {
	        throw new Error('Rule should declare lineBreaks: ' + regexp)
	      }

	      // store regex
	      parts.push(reCapture(pat));
	    }


	    // If there's no fallback rule, use the sticky flag so we only look for
	    // matches at the current index.
	    //
	    // If we don't support the sticky flag, then fake it using an irrefutable
	    // match (i.e. an empty pattern).
	    var fallbackRule = errorRule && errorRule.fallback;
	    var flags = hasSticky && !fallbackRule ? 'ym' : 'gm';
	    var suffix = hasSticky || fallbackRule ? '' : '|';

	    if (unicodeFlag === true) flags += "u";
	    var combined = new RegExp(reUnion(parts) + suffix, flags);
	    return {regexp: combined, groups: groups, fast: fast, error: errorRule || defaultErrorRule}
	  }

	  function compile(rules) {
	    var result = compileRules(toRules(rules));
	    return new Lexer({start: result}, 'start')
	  }

	  function checkStateGroup(g, name, map) {
	    var state = g && (g.push || g.next);
	    if (state && !map[state]) {
	      throw new Error("Missing state '" + state + "' (in token '" + g.defaultType + "' of state '" + name + "')")
	    }
	    if (g && g.pop && +g.pop !== 1) {
	      throw new Error("pop must be 1 (in token '" + g.defaultType + "' of state '" + name + "')")
	    }
	  }
	  function compileStates(states, start) {
	    var all = states.$all ? toRules(states.$all) : [];
	    delete states.$all;

	    var keys = Object.getOwnPropertyNames(states);
	    if (!start) start = keys[0];

	    var ruleMap = Object.create(null);
	    for (var i = 0; i < keys.length; i++) {
	      var key = keys[i];
	      ruleMap[key] = toRules(states[key]).concat(all);
	    }
	    for (var i = 0; i < keys.length; i++) {
	      var key = keys[i];
	      var rules = ruleMap[key];
	      var included = Object.create(null);
	      for (var j = 0; j < rules.length; j++) {
	        var rule = rules[j];
	        if (!rule.include) continue
	        var splice = [j, 1];
	        if (rule.include !== key && !included[rule.include]) {
	          included[rule.include] = true;
	          var newRules = ruleMap[rule.include];
	          if (!newRules) {
	            throw new Error("Cannot include nonexistent state '" + rule.include + "' (in state '" + key + "')")
	          }
	          for (var k = 0; k < newRules.length; k++) {
	            var newRule = newRules[k];
	            if (rules.indexOf(newRule) !== -1) continue
	            splice.push(newRule);
	          }
	        }
	        rules.splice.apply(rules, splice);
	        j--;
	      }
	    }

	    var map = Object.create(null);
	    for (var i = 0; i < keys.length; i++) {
	      var key = keys[i];
	      map[key] = compileRules(ruleMap[key], true);
	    }

	    for (var i = 0; i < keys.length; i++) {
	      var name = keys[i];
	      var state = map[name];
	      var groups = state.groups;
	      for (var j = 0; j < groups.length; j++) {
	        checkStateGroup(groups[j], name, map);
	      }
	      var fastKeys = Object.getOwnPropertyNames(state.fast);
	      for (var j = 0; j < fastKeys.length; j++) {
	        checkStateGroup(state.fast[fastKeys[j]], name, map);
	      }
	    }

	    return new Lexer(map, start)
	  }

	  function keywordTransform(map) {

	    // Use a JavaScript Map to map keywords to their corresponding token type
	    // unless Map is unsupported, then fall back to using an Object:
	    var isMap = typeof Map !== 'undefined';
	    var reverseMap = isMap ? new Map : Object.create(null);

	    var types = Object.getOwnPropertyNames(map);
	    for (var i = 0; i < types.length; i++) {
	      var tokenType = types[i];
	      var item = map[tokenType];
	      var keywordList = Array.isArray(item) ? item : [item];
	      keywordList.forEach(function(keyword) {
	        if (typeof keyword !== 'string') {
	          throw new Error("keyword must be string (in keyword '" + tokenType + "')")
	        }
	        if (isMap) {
	          reverseMap.set(keyword, tokenType);
	        } else {
	          reverseMap[keyword] = tokenType;
	        }
	      });
	    }
	    return function(k) {
	      return isMap ? reverseMap.get(k) : reverseMap[k]
	    }
	  }

	  /***************************************************************************/

	  var Lexer = function(states, state) {
	    this.startState = state;
	    this.states = states;
	    this.buffer = '';
	    this.stack = [];
	    this.reset();
	  };

	  Lexer.prototype.reset = function(data, info) {
	    this.buffer = data || '';
	    this.index = 0;
	    this.line = info ? info.line : 1;
	    this.col = info ? info.col : 1;
	    this.queuedToken = info ? info.queuedToken : null;
	    this.queuedText = info ? info.queuedText: "";
	    this.queuedThrow = info ? info.queuedThrow : null;
	    this.setState(info ? info.state : this.startState);
	    this.stack = info && info.stack ? info.stack.slice() : [];
	    return this
	  };

	  Lexer.prototype.save = function() {
	    return {
	      line: this.line,
	      col: this.col,
	      state: this.state,
	      stack: this.stack.slice(),
	      queuedToken: this.queuedToken,
	      queuedText: this.queuedText,
	      queuedThrow: this.queuedThrow,
	    }
	  };

	  Lexer.prototype.setState = function(state) {
	    if (!state || this.state === state) return
	    this.state = state;
	    var info = this.states[state];
	    this.groups = info.groups;
	    this.error = info.error;
	    this.re = info.regexp;
	    this.fast = info.fast;
	  };

	  Lexer.prototype.popState = function() {
	    this.setState(this.stack.pop());
	  };

	  Lexer.prototype.pushState = function(state) {
	    this.stack.push(this.state);
	    this.setState(state);
	  };

	  var eat = hasSticky ? function(re, buffer) { // assume re is /y
	    return re.exec(buffer)
	  } : function(re, buffer) { // assume re is /g
	    var match = re.exec(buffer);
	    // will always match, since we used the |(?:) trick
	    if (match[0].length === 0) {
	      return null
	    }
	    return match
	  };

	  Lexer.prototype._getGroup = function(match) {
	    var groupCount = this.groups.length;
	    for (var i = 0; i < groupCount; i++) {
	      if (match[i + 1] !== undefined) {
	        return this.groups[i]
	      }
	    }
	    throw new Error('Cannot find token type for matched text')
	  };

	  function tokenToString() {
	    return this.value
	  }

	  Lexer.prototype.next = function() {
	    var index = this.index;

	    // If a fallback token matched, we don't need to re-run the RegExp
	    if (this.queuedGroup) {
	      var token = this._token(this.queuedGroup, this.queuedText, index);
	      this.queuedGroup = null;
	      this.queuedText = "";
	      return token
	    }

	    var buffer = this.buffer;
	    if (index === buffer.length) {
	      return // EOF
	    }

	    // Fast matching for single characters
	    var group = this.fast[buffer.charCodeAt(index)];
	    if (group) {
	      return this._token(group, buffer.charAt(index), index)
	    }

	    // Execute RegExp
	    var re = this.re;
	    re.lastIndex = index;
	    var match = eat(re, buffer);

	    // Error tokens match the remaining buffer
	    var error = this.error;
	    if (match == null) {
	      return this._token(error, buffer.slice(index, buffer.length), index)
	    }

	    var group = this._getGroup(match);
	    var text = match[0];

	    if (error.fallback && match.index !== index) {
	      this.queuedGroup = group;
	      this.queuedText = text;

	      // Fallback tokens contain the unmatched portion of the buffer
	      return this._token(error, buffer.slice(index, match.index), index)
	    }

	    return this._token(group, text, index)
	  };

	  Lexer.prototype._token = function(group, text, offset) {
	    // count line breaks
	    var lineBreaks = 0;
	    if (group.lineBreaks) {
	      var matchNL = /\n/g;
	      var nl = 1;
	      if (text === '\n') {
	        lineBreaks = 1;
	      } else {
	        while (matchNL.exec(text)) { lineBreaks++; nl = matchNL.lastIndex; }
	      }
	    }

	    var token = {
	      type: (typeof group.type === 'function' && group.type(text)) || group.defaultType,
	      value: typeof group.value === 'function' ? group.value(text) : text,
	      text: text,
	      toString: tokenToString,
	      offset: offset,
	      lineBreaks: lineBreaks,
	      line: this.line,
	      col: this.col,
	    };
	    // nb. adding more props to token object will make V8 sad!

	    var size = text.length;
	    this.index += size;
	    this.line += lineBreaks;
	    if (lineBreaks !== 0) {
	      this.col = size - nl + 1;
	    } else {
	      this.col += size;
	    }

	    // throw, if no rule with {error: true}
	    if (group.shouldThrow) {
	      var err = new Error(this.formatError(token, "invalid syntax"));
	      throw err;
	    }

	    if (group.pop) this.popState();
	    else if (group.push) this.pushState(group.push);
	    else if (group.next) this.setState(group.next);

	    return token
	  };

	  if (typeof Symbol !== 'undefined' && Symbol.iterator) {
	    var LexerIterator = function(lexer) {
	      this.lexer = lexer;
	    };

	    LexerIterator.prototype.next = function() {
	      var token = this.lexer.next();
	      return {value: token, done: !token}
	    };

	    LexerIterator.prototype[Symbol.iterator] = function() {
	      return this
	    };

	    Lexer.prototype[Symbol.iterator] = function() {
	      return new LexerIterator(this)
	    };
	  }

	  Lexer.prototype.formatError = function(token, message) {
	    if (token == null) {
	      // An undefined token indicates EOF
	      var text = this.buffer.slice(this.index);
	      var token = {
	        text: text,
	        offset: this.index,
	        lineBreaks: text.indexOf('\n') === -1 ? 0 : 1,
	        line: this.line,
	        col: this.col,
	      };
	    }
	    
	    var numLinesAround = 2;
	    var firstDisplayedLine = Math.max(token.line - numLinesAround, 1);
	    var lastDisplayedLine = token.line + numLinesAround;
	    var lastLineDigits = String(lastDisplayedLine).length;
	    var displayedLines = lastNLines(
	        this.buffer, 
	        (this.line - token.line) + numLinesAround + 1
	      )
	      .slice(0, 5);
	    var errorLines = [];
	    errorLines.push(message + " at line " + token.line + " col " + token.col + ":");
	    errorLines.push("");
	    for (var i = 0; i < displayedLines.length; i++) {
	      var line = displayedLines[i];
	      var lineNo = firstDisplayedLine + i;
	      errorLines.push(pad(String(lineNo), lastLineDigits) + "  " + line);
	      if (lineNo === token.line) {
	        errorLines.push(pad("", lastLineDigits + token.col + 1) + "^");
	      }
	    }
	    return errorLines.join("\n")
	  };

	  Lexer.prototype.clone = function() {
	    return new Lexer(this.states, this.state)
	  };

	  Lexer.prototype.has = function(tokenType) {
	    return true
	  };


	  return {
	    compile: compile,
	    states: compileStates,
	    error: Object.freeze({error: true}),
	    fallback: Object.freeze({fallback: true}),
	    keywords: keywordTransform,
	  }

	})); 
} (moo$1));

var mooExports = moo$1.exports;

Object.defineProperty(parseley$1, '__esModule', { value: true });

var nearley = nearleyExports;
var moo = mooExports;

function _interopNamespace$1(e) {
    if (e && e.__esModule) return e;
    var n = Object.create(null);
    if (e) {
        Object.keys(e).forEach(function (k) {
            if (k !== 'default') {
                var d = Object.getOwnPropertyDescriptor(e, k);
                Object.defineProperty(n, k, d.get ? d : {
                    enumerable: true,
                    get: function () {
                        return e[k];
                    }
                });
            }
        });
    }
    n['default'] = e;
    return Object.freeze(n);
}

var moo__namespace = /*#__PURE__*/_interopNamespace$1(moo);

// Generated automatically by nearley, version 2.20.1
// http://github.com/Hardmath123/nearley
// Bypasses TS6133. Allow declared but unused functions.
// @ts-ignore
function id(d) { return d[0]; }
const lexer = moo__namespace.compile({
    ws: { match: /[ \t\r\n\f]+/, lineBreaks: true },
    idn: { match: /[a-zA-Z_-][a-zA-Z0-9_-]*/ },
    hashToken: { match: /#[a-zA-Z0-9_-]+/, value: (s) => s.slice(1) },
    str1: { match: /'(?:\\['\\]|[^\n'\\])*'/, value: (s) => s.slice(1, -1) },
    str2: { match: /"(?:\\["\\]|[^\n"\\])*"/, value: (s) => s.slice(1, -1) },
    asterisk: '*',
    fullstop: '.',
    comma: ',',
    lbr: '[',
    rbr: ']',
    eq: '=',
    gt: '>',
    vbar: '|',
    plus: '+',
    tilde: '~',
    caret: '^',
    dollar: '$',
    //colon:      ':',
    //lpar:       '(',
    //rpar:       ')',
});
function firstTokenValue(tokens) {
    return tokens[0].value;
}
function second(tokens) {
    return tokens[1];
}
function sumSpec([a0, a1, a2], [b0, b1, b2]) {
    return [a0 + b0, a1 + b1, a2 + b2];
}
const grammar = {
    Lexer: lexer,
    ParserRules: [
        { "name": "main", "symbols": ["_", "listSelector", "_"], "postprocess": second },
        { "name": "mainNoList", "symbols": ["_", "complexSelector", "_"], "postprocess": second },
        { "name": "listSelector", "symbols": ["complexSelector"], "postprocess": ([next]) => ({ type: 'list', list: [next] }) },
        { "name": "listSelector", "symbols": ["listSelector", "_", (lexer.has("comma") ? { type: "comma" } : comma), "_", "complexSelector"], "postprocess": ([acc, , , , next]) => ({ type: 'list', list: [...acc.list, next] }) },
        { "name": "complexSelector", "symbols": ["compoundSelector"], "postprocess": id },
        { "name": "complexSelector", "symbols": ["complexSelector", "__", "compoundSelector"], "postprocess": ([left, , right]) => ({
                type: 'compound',
                list: [...right.list, { type: 'combinator', combinator: ' ', left: left, specificity: left.specificity }],
                specificity: sumSpec(left.specificity, right.specificity)
            }) },
        { "name": "complexSelector", "symbols": ["complexSelector", "_", "combinator", "_", "compoundSelector"], "postprocess": ([left, , c, , right]) => ({
                type: 'compound',
                list: [...right.list, { type: 'combinator', combinator: c, left: left, specificity: left.specificity }],
                specificity: sumSpec(left.specificity, right.specificity)
            }) },
        { "name": "combinator", "symbols": [(lexer.has("gt") ? { type: "gt" } : gt)], "postprocess": () => '>' },
        { "name": "combinator", "symbols": [(lexer.has("plus") ? { type: "plus" } : plus)], "postprocess": () => '+' },
        { "name": "combinator", "symbols": [(lexer.has("tilde") ? { type: "tilde" } : tilde)], "postprocess": () => '~' },
        { "name": "combinator", "symbols": [(lexer.has("vbar") ? { type: "vbar" } : vbar), (lexer.has("vbar") ? { type: "vbar" } : vbar)], "postprocess": () => '||' },
        { "name": "compoundSelector", "symbols": ["typeSelector"], "postprocess": ([next]) => ({
                type: 'compound',
                list: [next],
                specificity: next.specificity
            }) },
        { "name": "compoundSelector", "symbols": ["subclassSelector"], "postprocess": ([next]) => ({
                type: 'compound',
                list: [next],
                specificity: next.specificity
            }) },
        { "name": "compoundSelector", "symbols": ["compoundSelector", "subclassSelector"], "postprocess": ([acc, next]) => ({
                type: 'compound',
                list: [...acc.list, next],
                specificity: sumSpec(acc.specificity, next.specificity)
            }) },
        { "name": "subclassSelector", "symbols": ["idSelector"], "postprocess": id },
        { "name": "subclassSelector", "symbols": ["classSelector"], "postprocess": id },
        { "name": "subclassSelector", "symbols": ["attrSelector"], "postprocess": id },
        { "name": "attrSelector", "symbols": ["attrPresenceSelector"], "postprocess": id },
        { "name": "attrSelector", "symbols": ["attrValueSelector"], "postprocess": id },
        { "name": "typeSelector", "symbols": ["tagSelector"], "postprocess": id },
        { "name": "typeSelector", "symbols": ["uniSelector"], "postprocess": id },
        { "name": "attrPresenceSelector", "symbols": [(lexer.has("lbr") ? { type: "lbr" } : lbr), "_", "wqname", "_", (lexer.has("rbr") ? { type: "rbr" } : rbr)], "postprocess": ([, , wqname]) => ({
                type: 'attrPresence',
                name: wqname.name,
                namespace: wqname.namespace,
                specificity: [0, 1, 0]
            })
        },
        { "name": "attrValueSelector", "symbols": [(lexer.has("lbr") ? { type: "lbr" } : lbr), "_", "wqname", "_", "attrMatcher", "_", "attrValue", "_", (lexer.has("rbr") ? { type: "rbr" } : rbr)], "postprocess": ([, , wqname, , matcher, , v]) => ({
                type: 'attrValue',
                name: wqname.name,
                namespace: wqname.namespace,
                matcher: matcher,
                value: v.value,
                modifier: v.modifier,
                specificity: [0, 1, 0]
            })
        },
        { "name": "attrMatcher", "symbols": [(lexer.has("eq") ? { type: "eq" } : eq)], "postprocess": () => '=' },
        { "name": "attrMatcher", "symbols": [(lexer.has("tilde") ? { type: "tilde" } : tilde), (lexer.has("eq") ? { type: "eq" } : eq)], "postprocess": () => '~=' },
        { "name": "attrMatcher", "symbols": [(lexer.has("vbar") ? { type: "vbar" } : vbar), (lexer.has("eq") ? { type: "eq" } : eq)], "postprocess": () => '|=' },
        { "name": "attrMatcher", "symbols": [(lexer.has("caret") ? { type: "caret" } : caret), (lexer.has("eq") ? { type: "eq" } : eq)], "postprocess": () => '^=' },
        { "name": "attrMatcher", "symbols": [(lexer.has("dollar") ? { type: "dollar" } : dollar), (lexer.has("eq") ? { type: "eq" } : eq)], "postprocess": () => '$=' },
        { "name": "attrMatcher", "symbols": [(lexer.has("asterisk") ? { type: "asterisk" } : asterisk), (lexer.has("eq") ? { type: "eq" } : eq)], "postprocess": () => '*=' },
        { "name": "attrValue", "symbols": ["str"], "postprocess": ([v]) => ({ value: v, modifier: null }) },
        { "name": "attrValue", "symbols": ["idn"], "postprocess": ([v]) => ({ value: v, modifier: null }) },
        { "name": "attrValue", "symbols": ["str", "_", "attrModifier"], "postprocess": ([v, , mod]) => ({ value: v, modifier: mod }) },
        { "name": "attrValue", "symbols": ["idn", "__", "attrModifier"], "postprocess": ([v, , mod]) => ({ value: v, modifier: mod }) },
        { "name": "attrModifier", "symbols": [{ "literal": "i" }], "postprocess": () => 'i' },
        { "name": "attrModifier", "symbols": [{ "literal": "I" }], "postprocess": () => 'i' },
        { "name": "attrModifier", "symbols": [{ "literal": "s" }], "postprocess": () => 's' },
        { "name": "attrModifier", "symbols": [{ "literal": "S" }], "postprocess": () => 's' },
        { "name": "idSelector", "symbols": [(lexer.has("hashToken") ? { type: "hashToken" } : hashToken)], "postprocess": ([{ value: name }]) => ({ type: 'id', name: name, specificity: [1, 0, 0] }) },
        { "name": "classSelector", "symbols": [(lexer.has("fullstop") ? { type: "fullstop" } : fullstop), "idn"], "postprocess": ([, name]) => ({ type: 'class', name: name, specificity: [0, 1, 0] }) },
        { "name": "tagSelector", "symbols": ["wqname"], "postprocess": ([wqname]) => ({
                type: 'tag',
                name: wqname.name,
                namespace: wqname.namespace,
                specificity: [0, 0, 1]
            })
        },
        { "name": "uniSelector", "symbols": [(lexer.has("asterisk") ? { type: "asterisk" } : asterisk)], "postprocess": () => ({ type: 'universal', namespace: null, specificity: [0, 0, 0] }) },
        { "name": "uniSelector", "symbols": ["ns", (lexer.has("asterisk") ? { type: "asterisk" } : asterisk)], "postprocess": ([ns]) => ({ type: 'universal', namespace: ns, specificity: [0, 0, 0] }) },
        { "name": "wqname", "symbols": ["idn"], "postprocess": ([name]) => ({ name: name, namespace: null }) },
        { "name": "wqname", "symbols": ["ns", "idn"], "postprocess": ([ns, name]) => ({ name: name, namespace: ns }) },
        { "name": "ns", "symbols": [(lexer.has("vbar") ? { type: "vbar" } : vbar)], "postprocess": () => '' },
        { "name": "ns", "symbols": ["idn", (lexer.has("vbar") ? { type: "vbar" } : vbar)], "postprocess": id },
        { "name": "str", "symbols": [(lexer.has("str1") ? { type: "str1" } : str1)], "postprocess": firstTokenValue },
        { "name": "str", "symbols": [(lexer.has("str2") ? { type: "str2" } : str2)], "postprocess": firstTokenValue },
        { "name": "idn", "symbols": [(lexer.has("idn") ? { type: "idn" } : idn)], "postprocess": firstTokenValue },
        { "name": "_$ebnf$1", "symbols": [(lexer.has("ws") ? { type: "ws" } : ws)], "postprocess": id },
        { "name": "_$ebnf$1", "symbols": [], "postprocess": () => null },
        { "name": "_", "symbols": ["_$ebnf$1"], "postprocess": () => null },
        { "name": "__", "symbols": [(lexer.has("ws") ? { type: "ws" } : ws)], "postprocess": () => null }
    ],
    ParserStart: "main",
};

var ast$2 = /*#__PURE__*/Object.freeze({
    __proto__: null
});

// Passing the start argument to a parser or grammar constructor
// doesn't seem to work as expected.
const compiledRulesNoList = { ...grammar, ParserStart: 'mainNoList' };
/**
 * Parse a CSS selector string.
 *
 * This function supports comma-separated selector lists
 * and always returns an AST starting from a node of type `list`.
 *
 * @param str - CSS selector string (can contain commas).
 */
function parse(str) {
    return _parse(grammar, str);
}
/**
 * Parse a CSS selector string.
 *
 * This function does not support comma-separated selector lists
 * and always returns an AST starting from a node of type `compound`.
 *
 * @param str - CSS selector string (no commas).
 */
function parse1(str) {
    return _parse(compiledRulesNoList, str);
}
function _parse(compiledRules1, str) {
    const parser = new nearley.Parser(nearley.Grammar.fromCompiled(compiledRules1));
    parser.feed(str);
    if (parser.results.length === 0) {
        throw new Error('Failed to parse - input string might be incomplete.');
    }
    return parser.results[0];
}
/**
 * Convert a selector AST back to a string representation.
 *
 * Note: formatting is not preserved in the AST.
 *
 * @param selector - A selector AST object.
 */
function serialize(selector) {
    if (!selector.type) {
        throw new Error('This is not an AST node.');
    }
    switch (selector.type) {
        case 'universal':
            return _serNs(selector.namespace) + '*';
        case 'tag':
            return _serNs(selector.namespace) + selector.name;
        case 'class':
            return '.' + selector.name;
        case 'id':
            return '#' + selector.name;
        case 'attrPresence':
            return `[${_serNs(selector.namespace)}${selector.name}]`;
        case 'attrValue':
            return `[${_serNs(selector.namespace)}${selector.name}${selector.matcher}${_serStr(selector.value)}${(selector.modifier ? selector.modifier : '')}]`;
        case 'combinator':
            return serialize(selector.left) + selector.combinator;
        case 'compound':
            return selector.list.reduce((acc, node) => {
                if (node.type === 'combinator') {
                    return serialize(node) + acc;
                }
                else {
                    return acc + serialize(node);
                }
            }, '');
        case 'list':
            return selector.list.map(serialize).join(',');
    }
}
function _serNs(ns) {
    return (ns || ns === '')
        ? ns + '|'
        : '';
}
function _serStr(str) {
    if (str.indexOf('"') === -1) {
        return `"${str}"`;
    }
    else if (str.indexOf("'") === -1) {
        return `'${str}'`;
    }
    else {
        return `"${str.replace('"', '\\"')}"`;
    }
}
/**
 * Modifies the given AST **in place** to have all internal arrays
 * in a stable order. Returns the AST.
 *
 * Intended for consitent processing and normalized `serialize()` output.
 *
 * @param selector - A selector AST object.
 */
function normalize(selector) {
    if (!selector.type) {
        throw new Error('This is not an AST node.');
    }
    switch (selector.type) {
        case 'compound': {
            selector.list.forEach(normalize);
            selector.list.sort((a, b) => _compareArrays(_getSelectorPriority(a), _getSelectorPriority(b)));
            break;
        }
        case 'combinator': {
            normalize(selector.left);
            break;
        }
        case 'list': {
            selector.list.forEach(normalize);
            selector.list.sort((a, b) => (serialize(a) < serialize(b)) ? -1 : 1);
            break;
        }
    }
    return selector;
}
function _getSelectorPriority(selector) {
    switch (selector.type) {
        case 'universal':
            return [1];
        case 'tag':
            return [1];
        case 'id':
            return [2];
        case 'class':
            return [3, selector.name];
        case 'attrPresence':
            return [4, serialize(selector)];
        case 'attrValue':
            return [5, serialize(selector)];
        case 'combinator':
            return [15, serialize(selector)];
    }
}
/**
 * Compare selectors based on their specificity.
 *
 * Usable as a comparator for sorting.
 *
 * @param a - First selector.
 * @param b - Second selector.
 */
function compareSelectors(a, b) {
    return _compareArrays(a.specificity, b.specificity);
}
/**
 * Compare specificity values without reducing them
 * as arbitrary base numbers.
 *
 * Usable as a comparator for sorting.
 *
 * @param a - First specificity value.
 * @param b - Second specificity value.
 */
function compareSpecificity(a, b) {
    return _compareArrays(a, b);
}
function _compareArrays(a, b) {
    if (!Array.isArray(a) || !Array.isArray(b)) {
        throw new Error('Arguments must be arrays.');
    }
    const shorter = (a.length < b.length) ? a.length : b.length;
    for (let i = 0; i < shorter; i++) {
        if (a[i] === b[i]) {
            continue;
        }
        return (a[i] < b[i]) ? -1 : 1;
    }
    return a.length - b.length;
}

parseley$1.Ast = ast$2;
parseley$1.compareSelectors = compareSelectors;
parseley$1.compareSpecificity = compareSpecificity;
parseley$1.normalize = normalize;
parseley$1.parse = parse;
parseley$1.parse1 = parse1;
parseley$1.serialize = serialize;

Object.defineProperty(selderee$2, '__esModule', { value: true });

var parseley = parseley$1;

function _interopNamespace(e) {
    if (e && e.__esModule) return e;
    var n = Object.create(null);
    if (e) {
        Object.keys(e).forEach(function (k) {
            if (k !== 'default') {
                var d = Object.getOwnPropertyDescriptor(e, k);
                Object.defineProperty(n, k, d.get ? d : {
                    enumerable: true,
                    get: function () {
                        return e[k];
                    }
                });
            }
        });
    }
    n['default'] = e;
    return Object.freeze(n);
}

var parseley__namespace = /*#__PURE__*/_interopNamespace(parseley);

var Ast = /*#__PURE__*/Object.freeze({
    __proto__: null
});

var Types = /*#__PURE__*/Object.freeze({
    __proto__: null
});

/**
 * A {@link BuilderFunction} implementation.
 *
 * Produces a string representation of the tree
 * for testing and debug purposes.
 *
 * Only accepts `string` as the associated value type.
 * Map your input collection before creating a {@link DecisionTree}
 * if you want to use it with a different type -
 * the decision on how to stringify the value is up to you.
 *
 * @param nodes - nodes from the root level of the decision tree.
 * @returns the string representation of the tree.
 */
const treeify = (nodes) => '▽\n' + treeifyArray(nodes, thinLines);
const thinLines = [['├─', '│ '], ['└─', '  ']];
const heavyLines = [['┠─', '┃ '], ['┖─', '  ']];
const doubleLines = [['╟─', '║ '], ['╙─', '  ']];
function treeifyArray(nodes, tpl = heavyLines) {
    return prefixItems(tpl, nodes.map(n => treeifyNode(n)));
}
function treeifyNode(node) {
    switch (node.type) {
        case 'terminal': {
            const vctr = node.valueContainer;
            return `◁ #${vctr.index} ${JSON.stringify(vctr.specificity)} ${vctr.value}`;
        }
        case 'tagName':
            return `◻ Tag name\n${treeifyArray(node.variants, doubleLines)}`;
        case 'attrValue':
            return `▣ Attr value: ${node.name}\n${treeifyArray(node.matchers, doubleLines)}`;
        case 'attrPresence':
            return `◨ Attr presence: ${node.name}\n${treeifyArray(node.cont)}`;
        case 'pushElement':
            return `◉ Push element: ${node.combinator}\n${treeifyArray(node.cont, thinLines)}`;
        case 'popElement':
            return `◌ Pop element\n${treeifyArray(node.cont, thinLines)}`;
        case 'variant':
            return `◇ = ${node.value}\n${treeifyArray(node.cont)}`;
        case 'matcher':
            return `◈ ${node.matcher} "${node.value}"${node.modifier || ''}\n${treeifyArray(node.cont)}`;
    }
}
function prefixItems(tpl, items) {
    return items
        .map((item, i, { length }) => prefixItem(tpl, item, i === length - 1))
        .join('\n');
}
function prefixItem(tpl, item, tail = true) {
    const tpl1 = tpl[tail ? 1 : 0];
    return tpl1[0] + item.split('\n').join('\n' + tpl1[1]);
}

var TreeifyBuilder = /*#__PURE__*/Object.freeze({
    __proto__: null,
    treeify: treeify
});

/**
 * CSS selectors decision tree.
 * Data structure that weaves similar selectors together
 * in order to minimize the number of checks required
 * to find the ones matching a given HTML element.
 *
 * Converted into a functioning implementation via plugins
 * tailored for specific DOM ASTs.
 *
 * @typeParam V - the type of values associated with selectors.
 */
class DecisionTree {
    /**
     * Create new DecisionTree object.
     *
     * @param input - an array containing all selectors
     * paired with associated values.
     *
     * @typeParam V - the type of values associated with selectors.
     */
    constructor(input) {
        this.branches = weave(toAstTerminalPairs(input));
    }
    /**
     * Turn this decision tree into a usable form.
     *
     * @typeParam R - return type defined by the builder function.
     *
     * @param builder - the builder function.
     *
     * @returns the decision tree in a form ready for use.
     */
    build(builder) {
        return builder(this.branches);
    }
}
function toAstTerminalPairs(array) {
    const len = array.length;
    const results = new Array(len);
    for (let i = 0; i < len; i++) {
        const [selectorString, val] = array[i];
        const ast = preprocess(parseley__namespace.parse1(selectorString));
        results[i] = {
            ast: ast,
            terminal: {
                type: 'terminal',
                valueContainer: { index: i, value: val, specificity: ast.specificity }
            }
        };
    }
    return results;
}
function preprocess(ast) {
    reduceSelectorVariants(ast);
    parseley__namespace.normalize(ast);
    return ast;
}
function reduceSelectorVariants(ast) {
    const newList = [];
    ast.list.forEach(sel => {
        switch (sel.type) {
            case 'class':
                newList.push({
                    matcher: '~=',
                    modifier: null,
                    name: 'class',
                    namespace: null,
                    specificity: sel.specificity,
                    type: 'attrValue',
                    value: sel.name,
                });
                break;
            case 'id':
                newList.push({
                    matcher: '=',
                    modifier: null,
                    name: 'id',
                    namespace: null,
                    specificity: sel.specificity,
                    type: 'attrValue',
                    value: sel.name,
                });
                break;
            case 'combinator':
                reduceSelectorVariants(sel.left);
                newList.push(sel);
                break;
            case 'universal':
                // skip it
                break;
            default:
                newList.push(sel);
                break;
        }
    });
    ast.list = newList;
}
function weave(items) {
    const branches = [];
    while (items.length) {
        const topKind = findTopKey(items, (sel) => true, getSelectorKind);
        const { matches, nonmatches, empty } = breakByKind(items, topKind);
        items = nonmatches;
        if (matches.length) {
            branches.push(branchOfKind(topKind, matches));
        }
        if (empty.length) {
            branches.push(...terminate(empty));
        }
    }
    return branches;
}
function terminate(items) {
    const results = [];
    for (const item of items) {
        const terminal = item.terminal;
        if (terminal.type === 'terminal') {
            results.push(terminal);
        }
        else { // popElement - lift contained terminals
            const { matches, rest } = partition(terminal.cont, (node) => node.type === 'terminal');
            matches.forEach((node) => results.push(node));
            if (rest.length) {
                terminal.cont = rest;
                results.push(terminal);
            }
        }
    }
    return results;
}
function breakByKind(items, selectedKind) {
    const matches = [];
    const nonmatches = [];
    const empty = [];
    for (const item of items) {
        const simpsels = item.ast.list;
        if (simpsels.length) {
            const isMatch = simpsels.some(node => getSelectorKind(node) === selectedKind);
            (isMatch ? matches : nonmatches).push(item);
        }
        else {
            empty.push(item);
        }
    }
    return { matches, nonmatches, empty };
}
function getSelectorKind(sel) {
    switch (sel.type) {
        case 'attrPresence':
            return `attrPresence ${sel.name}`;
        case 'attrValue':
            return `attrValue ${sel.name}`;
        case 'combinator':
            return `combinator ${sel.combinator}`;
        default:
            return sel.type;
    }
}
function branchOfKind(kind, items) {
    if (kind === 'tag') {
        return tagNameBranch(items);
    }
    if (kind.startsWith('attrValue ')) {
        return attrValueBranch(kind.substring(10), items);
    }
    if (kind.startsWith('attrPresence ')) {
        return attrPresenceBranch(kind.substring(13), items);
    }
    if (kind === 'combinator >') {
        return combinatorBranch('>', items);
    }
    if (kind === 'combinator +') {
        return combinatorBranch('+', items);
    }
    throw new Error(`Unsupported selector kind: ${kind}`);
}
function tagNameBranch(items) {
    const groups = spliceAndGroup(items, (x) => x.type === 'tag', (x) => x.name);
    const variants = Object.entries(groups).map(([name, group]) => ({
        type: 'variant',
        value: name,
        cont: weave(group.items)
    }));
    return {
        type: 'tagName',
        variants: variants
    };
}
function attrPresenceBranch(name, items) {
    for (const item of items) {
        spliceSimpleSelector(item, (x) => (x.type === 'attrPresence') && (x.name === name));
    }
    return {
        type: 'attrPresence',
        name: name,
        cont: weave(items)
    };
}
function attrValueBranch(name, items) {
    const groups = spliceAndGroup(items, (x) => (x.type === 'attrValue') && (x.name === name), (x) => `${x.matcher} ${x.modifier || ''} ${x.value}`);
    const matchers = [];
    for (const group of Object.values(groups)) {
        const sel = group.oneSimpleSelector;
        const predicate = getAttrPredicate(sel);
        const continuation = weave(group.items);
        matchers.push({
            type: 'matcher',
            matcher: sel.matcher,
            modifier: sel.modifier,
            value: sel.value,
            predicate: predicate,
            cont: continuation
        });
    }
    return {
        type: 'attrValue',
        name: name,
        matchers: matchers
    };
}
function getAttrPredicate(sel) {
    if (sel.modifier === 'i') {
        const expected = sel.value.toLowerCase();
        switch (sel.matcher) {
            case '=':
                return (actual) => expected === actual.toLowerCase();
            case '~=':
                return (actual) => actual.toLowerCase().split(/[ \t]+/).includes(expected);
            case '^=':
                return (actual) => actual.toLowerCase().startsWith(expected);
            case '$=':
                return (actual) => actual.toLowerCase().endsWith(expected);
            case '*=':
                return (actual) => actual.toLowerCase().includes(expected);
            case '|=':
                return (actual) => {
                    const lower = actual.toLowerCase();
                    return (expected === lower) || (lower.startsWith(expected) && lower[expected.length] === '-');
                };
        }
    }
    else {
        const expected = sel.value;
        switch (sel.matcher) {
            case '=':
                return (actual) => expected === actual;
            case '~=':
                return (actual) => actual.split(/[ \t]+/).includes(expected);
            case '^=':
                return (actual) => actual.startsWith(expected);
            case '$=':
                return (actual) => actual.endsWith(expected);
            case '*=':
                return (actual) => actual.includes(expected);
            case '|=':
                return (actual) => (expected === actual) || (actual.startsWith(expected) && actual[expected.length] === '-');
        }
    }
}
function combinatorBranch(combinator, items) {
    const groups = spliceAndGroup(items, (x) => (x.type === 'combinator') && (x.combinator === combinator), (x) => parseley__namespace.serialize(x.left));
    const leftItems = [];
    for (const group of Object.values(groups)) {
        const rightCont = weave(group.items);
        const leftAst = group.oneSimpleSelector.left;
        leftItems.push({
            ast: leftAst,
            terminal: { type: 'popElement', cont: rightCont }
        });
    }
    return {
        type: 'pushElement',
        combinator: combinator,
        cont: weave(leftItems)
    };
}
function spliceAndGroup(items, predicate, keyCallback) {
    const groups = {};
    while (items.length) {
        const bestKey = findTopKey(items, predicate, keyCallback);
        const bestKeyPredicate = (sel) => predicate(sel) && keyCallback(sel) === bestKey;
        const hasBestKeyPredicate = (item) => item.ast.list.some(bestKeyPredicate);
        const { matches, rest } = partition1(items, hasBestKeyPredicate);
        let oneSimpleSelector = null;
        for (const item of matches) {
            const splicedNode = spliceSimpleSelector(item, bestKeyPredicate);
            if (!oneSimpleSelector) {
                oneSimpleSelector = splicedNode;
            }
        }
        if (oneSimpleSelector == null) {
            throw new Error('No simple selector is found.');
        }
        groups[bestKey] = { oneSimpleSelector: oneSimpleSelector, items: matches };
        items = rest;
    }
    return groups;
}
function spliceSimpleSelector(item, predicate) {
    const simpsels = item.ast.list;
    const matches = new Array(simpsels.length);
    let firstIndex = -1;
    for (let i = simpsels.length; i-- > 0;) {
        if (predicate(simpsels[i])) {
            matches[i] = true;
            firstIndex = i;
        }
    }
    if (firstIndex == -1) {
        throw new Error(`Couldn't find the required simple selector.`);
    }
    const result = simpsels[firstIndex];
    item.ast.list = simpsels.filter((sel, i) => !matches[i]);
    return result;
}
function findTopKey(items, predicate, keyCallback) {
    const candidates = {};
    for (const item of items) {
        const candidates1 = {};
        for (const node of item.ast.list.filter(predicate)) {
            candidates1[keyCallback(node)] = true;
        }
        for (const key of Object.keys(candidates1)) {
            if (candidates[key]) {
                candidates[key]++;
            }
            else {
                candidates[key] = 1;
            }
        }
    }
    let topKind = '';
    let topCounter = 0;
    for (const entry of Object.entries(candidates)) {
        if (entry[1] > topCounter) {
            topKind = entry[0];
            topCounter = entry[1];
        }
    }
    return topKind;
}
function partition(src, predicate) {
    const matches = [];
    const rest = [];
    for (const x of src) {
        if (predicate(x)) {
            matches.push(x);
        }
        else {
            rest.push(x);
        }
    }
    return { matches, rest };
}
function partition1(src, predicate) {
    const matches = [];
    const rest = [];
    for (const x of src) {
        if (predicate(x)) {
            matches.push(x);
        }
        else {
            rest.push(x);
        }
    }
    return { matches, rest };
}

/**
 * Simple wrapper around the matcher function.
 * Recommended return type for builder plugins.
 *
 * @typeParam L - the type of HTML Element in the targeted DOM AST.
 * @typeParam V - the type of associated values.
 */
class Picker {
    /**
     * Create new Picker object.
     *
     * @typeParam L - the type of HTML Element in the targeted DOM AST.
     * @typeParam V - the type of associated values.
     *
     * @param f - the function that matches an element
     * and returns all associated values.
     */
    constructor(f) {
        this.f = f;
    }
    /**
     * Run the selectors decision tree against one HTML Element
     * and return all matched associated values
     * along with selector specificities.
     *
     * Client code then decides how to further process them
     * (sort, filter, etc).
     *
     * @param el - an HTML Element.
     *
     * @returns all associated values along with
     * selector specificities for all matched selectors.
     */
    pickAll(el) {
        return this.f(el);
    }
    /**
     * Run the selectors decision tree against one HTML Element
     * and choose the value from the most specific mached selector.
     *
     * @param el - an HTML Element.
     *
     * @param preferFirst - option to define which value to choose
     * when there are multiple matches with equal specificity.
     *
     * @returns the value from the most specific mached selector
     * or `null` if nothing matched.
     */
    pick1(el, preferFirst = false) {
        const results = this.f(el);
        const len = results.length;
        if (len === 0) {
            return null;
        }
        if (len === 1) {
            return results[0].value;
        }
        const comparator = (preferFirst)
            ? comparatorPreferFirst
            : comparatorPreferLast;
        let result = results[0];
        for (let i = 1; i < len; i++) {
            const next = results[i];
            if (comparator(result, next)) {
                result = next;
            }
        }
        return result.value;
    }
}
function comparatorPreferFirst(acc, next) {
    const diff = parseley.compareSpecificity(next.specificity, acc.specificity);
    return diff > 0 || (diff === 0 && next.index < acc.index);
}
function comparatorPreferLast(acc, next) {
    const diff = parseley.compareSpecificity(next.specificity, acc.specificity);
    return diff > 0 || (diff === 0 && next.index > acc.index);
}

selderee$2.Ast = Ast;
selderee$2.DecisionTree = DecisionTree;
selderee$2.Picker = Picker;
selderee$2.Treeify = TreeifyBuilder;
selderee$2.Types = Types;

Object.defineProperty(hp2Builder$2, '__esModule', { value: true });

var domhandler = lib$5;
var selderee$1 = selderee$2;

/**
 * A {@link BuilderFunction} implementation.
 *
 * Creates a function (in a {@link Picker} wrapper) that can run
 * the decision tree against `htmlparser2` `Element` nodes.
 *
 * @typeParam V - the type of values associated with selectors.
 *
 * @param nodes - nodes ({@link DecisionTreeNode})
 * from the root level of the decision tree.
 *
 * @returns a {@link Picker} object.
 */
function hp2Builder$1(nodes) {
    return new selderee$1.Picker(handleArray(nodes));
}
// ==============================================
function handleArray(nodes) {
    const matchers = nodes.map(handleNode);
    return (el, ...tail) => flatMap(matchers, m => m(el, ...tail));
}
function handleNode(node) {
    switch (node.type) {
        case 'terminal': {
            const result = [node.valueContainer];
            return (el, ...tail) => result;
        }
        case 'tagName':
            return handleTagName(node);
        case 'attrValue':
            return handleAttrValueName(node);
        case 'attrPresence':
            return handleAttrPresenceName(node);
        case 'pushElement':
            return handlePushElementNode(node);
        case 'popElement':
            return handlePopElementNode(node);
    }
}
function handleTagName(node) {
    const variants = {};
    for (const variant of node.variants) {
        variants[variant.value] = handleArray(variant.cont);
    }
    return (el, ...tail) => {
        const continuation = variants[el.name];
        return (continuation) ? continuation(el, ...tail) : [];
    };
}
function handleAttrPresenceName(node) {
    const attrName = node.name;
    const continuation = handleArray(node.cont);
    return (el, ...tail) => (Object.prototype.hasOwnProperty.call(el.attribs, attrName))
        ? continuation(el, ...tail)
        : [];
}
function handleAttrValueName(node) {
    const callbacks = [];
    for (const matcher of node.matchers) {
        const predicate = matcher.predicate;
        const continuation = handleArray(matcher.cont);
        callbacks.push((attr, el, ...tail) => (predicate(attr) ? continuation(el, ...tail) : []));
    }
    const attrName = node.name;
    return (el, ...tail) => {
        const attr = el.attribs[attrName];
        return (attr || attr === '')
            ? flatMap(callbacks, cb => cb(attr, el, ...tail))
            : [];
    };
}
function handlePushElementNode(node) {
    const continuation = handleArray(node.cont);
    const leftElementGetter = (node.combinator === '+')
        ? getPrecedingElement
        : getParentElement;
    return (el, ...tail) => {
        const next = leftElementGetter(el);
        if (next === null) {
            return [];
        }
        return continuation(next, el, ...tail);
    };
}
const getPrecedingElement = (el) => {
    const prev = el.prev;
    if (prev === null) {
        return null;
    }
    return (domhandler.isTag(prev)) ? prev : getPrecedingElement(prev);
};
const getParentElement = (el) => {
    const parent = el.parent;
    return (parent && domhandler.isTag(parent)) ? parent : null;
};
function handlePopElementNode(node) {
    const continuation = handleArray(node.cont);
    return (el, next, ...tail) => continuation(next, ...tail);
}
// Can be removed after transition to Node 12.
function flatMap(items, mapper) {
    return [].concat(...amap(items, mapper));
}
function amap(items, mapper) {
    const len = items.length;
    const res = new Array(len);
    for (let i = 0; i < len; i++) {
        res[i] = mapper(items[i]);
    }
    return res;
}

hp2Builder$2.hp2Builder = hp2Builder$1;

var isMergeableObject = function isMergeableObject(value) {
	return isNonNullObject(value)
		&& !isSpecial(value)
};

function isNonNullObject(value) {
	return !!value && typeof value === 'object'
}

function isSpecial(value) {
	var stringValue = Object.prototype.toString.call(value);

	return stringValue === '[object RegExp]'
		|| stringValue === '[object Date]'
		|| isReactElement(value)
}

// see https://github.com/facebook/react/blob/b5ac963fb791d1298e7f396236383bc955f916c1/src/isomorphic/classic/element/ReactElement.js#L21-L25
var canUseSymbol = typeof Symbol === 'function' && Symbol.for;
var REACT_ELEMENT_TYPE = canUseSymbol ? Symbol.for('react.element') : 0xeac7;

function isReactElement(value) {
	return value.$$typeof === REACT_ELEMENT_TYPE
}

function emptyTarget(val) {
	return Array.isArray(val) ? [] : {}
}

function cloneUnlessOtherwiseSpecified(value, options) {
	return (options.clone !== false && options.isMergeableObject(value))
		? deepmerge(emptyTarget(value), value, options)
		: value
}

function defaultArrayMerge(target, source, options) {
	return target.concat(source).map(function(element) {
		return cloneUnlessOtherwiseSpecified(element, options)
	})
}

function getMergeFunction(key, options) {
	if (!options.customMerge) {
		return deepmerge
	}
	var customMerge = options.customMerge(key);
	return typeof customMerge === 'function' ? customMerge : deepmerge
}

function getEnumerableOwnPropertySymbols(target) {
	return Object.getOwnPropertySymbols
		? Object.getOwnPropertySymbols(target).filter(function(symbol) {
			return Object.propertyIsEnumerable.call(target, symbol)
		})
		: []
}

function getKeys(target) {
	return Object.keys(target).concat(getEnumerableOwnPropertySymbols(target))
}

function propertyIsOnObject(object, property) {
	try {
		return property in object
	} catch(_) {
		return false
	}
}

// Protects from prototype poisoning and unexpected merging up the prototype chain.
function propertyIsUnsafe(target, key) {
	return propertyIsOnObject(target, key) // Properties are safe to merge if they don't exist in the target yet,
		&& !(Object.hasOwnProperty.call(target, key) // unsafe if they exist up the prototype chain,
			&& Object.propertyIsEnumerable.call(target, key)) // and also unsafe if they're nonenumerable.
}

function mergeObject(target, source, options) {
	var destination = {};
	if (options.isMergeableObject(target)) {
		getKeys(target).forEach(function(key) {
			destination[key] = cloneUnlessOtherwiseSpecified(target[key], options);
		});
	}
	getKeys(source).forEach(function(key) {
		if (propertyIsUnsafe(target, key)) {
			return
		}

		if (propertyIsOnObject(target, key) && options.isMergeableObject(source[key])) {
			destination[key] = getMergeFunction(key, options)(target[key], source[key], options);
		} else {
			destination[key] = cloneUnlessOtherwiseSpecified(source[key], options);
		}
	});
	return destination
}

function deepmerge(target, source, options) {
	options = options || {};
	options.arrayMerge = options.arrayMerge || defaultArrayMerge;
	options.isMergeableObject = options.isMergeableObject || isMergeableObject;
	// cloneUnlessOtherwiseSpecified is added to `options` so that custom arrayMerge()
	// implementations can use it. The caller may not replace it.
	options.cloneUnlessOtherwiseSpecified = cloneUnlessOtherwiseSpecified;

	var sourceIsArray = Array.isArray(source);
	var targetIsArray = Array.isArray(target);
	var sourceAndTargetTypesMatch = sourceIsArray === targetIsArray;

	if (!sourceAndTargetTypesMatch) {
		return cloneUnlessOtherwiseSpecified(source, options)
	} else if (sourceIsArray) {
		return options.arrayMerge(target, source, options)
	} else {
		return mergeObject(target, source, options)
	}
}

deepmerge.all = function deepmergeAll(array, options) {
	if (!Array.isArray(array)) {
		throw new Error('first argument should be an array')
	}

	return array.reduce(function(prev, next) {
		return deepmerge(prev, next, options)
	}, {})
};

var deepmerge_1 = deepmerge;

var cjs = deepmerge_1;

var he$2 = {exports: {}};

/*! https://mths.be/he v1.2.0 by @mathias | MIT license */
he$2.exports;

(function (module, exports) {
(function(root) {

		// Detect free variables `exports`.
		var freeExports = exports;

		// Detect free variable `module`.
		var freeModule = module &&
			module.exports == freeExports && module;

		// Detect free variable `global`, from Node.js or Browserified code,
		// and use it as `root`.
		var freeGlobal = typeof commonjsGlobal == 'object' && commonjsGlobal;
		if (freeGlobal.global === freeGlobal || freeGlobal.window === freeGlobal) {
			root = freeGlobal;
		}

		/*--------------------------------------------------------------------------*/

		// All astral symbols.
		var regexAstralSymbols = /[\uD800-\uDBFF][\uDC00-\uDFFF]/g;
		// All ASCII symbols (not just printable ASCII) except those listed in the
		// first column of the overrides table.
		// https://html.spec.whatwg.org/multipage/syntax.html#table-charref-overrides
		var regexAsciiWhitelist = /[\x01-\x7F]/g;
		// All BMP symbols that are not ASCII newlines, printable ASCII symbols, or
		// code points listed in the first column of the overrides table on
		// https://html.spec.whatwg.org/multipage/syntax.html#table-charref-overrides.
		var regexBmpWhitelist = /[\x01-\t\x0B\f\x0E-\x1F\x7F\x81\x8D\x8F\x90\x9D\xA0-\uFFFF]/g;

		var regexEncodeNonAscii = /<\u20D2|=\u20E5|>\u20D2|\u205F\u200A|\u219D\u0338|\u2202\u0338|\u2220\u20D2|\u2229\uFE00|\u222A\uFE00|\u223C\u20D2|\u223D\u0331|\u223E\u0333|\u2242\u0338|\u224B\u0338|\u224D\u20D2|\u224E\u0338|\u224F\u0338|\u2250\u0338|\u2261\u20E5|\u2264\u20D2|\u2265\u20D2|\u2266\u0338|\u2267\u0338|\u2268\uFE00|\u2269\uFE00|\u226A\u0338|\u226A\u20D2|\u226B\u0338|\u226B\u20D2|\u227F\u0338|\u2282\u20D2|\u2283\u20D2|\u228A\uFE00|\u228B\uFE00|\u228F\u0338|\u2290\u0338|\u2293\uFE00|\u2294\uFE00|\u22B4\u20D2|\u22B5\u20D2|\u22D8\u0338|\u22D9\u0338|\u22DA\uFE00|\u22DB\uFE00|\u22F5\u0338|\u22F9\u0338|\u2933\u0338|\u29CF\u0338|\u29D0\u0338|\u2A6D\u0338|\u2A70\u0338|\u2A7D\u0338|\u2A7E\u0338|\u2AA1\u0338|\u2AA2\u0338|\u2AAC\uFE00|\u2AAD\uFE00|\u2AAF\u0338|\u2AB0\u0338|\u2AC5\u0338|\u2AC6\u0338|\u2ACB\uFE00|\u2ACC\uFE00|\u2AFD\u20E5|[\xA0-\u0113\u0116-\u0122\u0124-\u012B\u012E-\u014D\u0150-\u017E\u0192\u01B5\u01F5\u0237\u02C6\u02C7\u02D8-\u02DD\u0311\u0391-\u03A1\u03A3-\u03A9\u03B1-\u03C9\u03D1\u03D2\u03D5\u03D6\u03DC\u03DD\u03F0\u03F1\u03F5\u03F6\u0401-\u040C\u040E-\u044F\u0451-\u045C\u045E\u045F\u2002-\u2005\u2007-\u2010\u2013-\u2016\u2018-\u201A\u201C-\u201E\u2020-\u2022\u2025\u2026\u2030-\u2035\u2039\u203A\u203E\u2041\u2043\u2044\u204F\u2057\u205F-\u2063\u20AC\u20DB\u20DC\u2102\u2105\u210A-\u2113\u2115-\u211E\u2122\u2124\u2127-\u2129\u212C\u212D\u212F-\u2131\u2133-\u2138\u2145-\u2148\u2153-\u215E\u2190-\u219B\u219D-\u21A7\u21A9-\u21AE\u21B0-\u21B3\u21B5-\u21B7\u21BA-\u21DB\u21DD\u21E4\u21E5\u21F5\u21FD-\u2205\u2207-\u2209\u220B\u220C\u220F-\u2214\u2216-\u2218\u221A\u221D-\u2238\u223A-\u2257\u2259\u225A\u225C\u225F-\u2262\u2264-\u228B\u228D-\u229B\u229D-\u22A5\u22A7-\u22B0\u22B2-\u22BB\u22BD-\u22DB\u22DE-\u22E3\u22E6-\u22F7\u22F9-\u22FE\u2305\u2306\u2308-\u2310\u2312\u2313\u2315\u2316\u231C-\u231F\u2322\u2323\u232D\u232E\u2336\u233D\u233F\u237C\u23B0\u23B1\u23B4-\u23B6\u23DC-\u23DF\u23E2\u23E7\u2423\u24C8\u2500\u2502\u250C\u2510\u2514\u2518\u251C\u2524\u252C\u2534\u253C\u2550-\u256C\u2580\u2584\u2588\u2591-\u2593\u25A1\u25AA\u25AB\u25AD\u25AE\u25B1\u25B3-\u25B5\u25B8\u25B9\u25BD-\u25BF\u25C2\u25C3\u25CA\u25CB\u25EC\u25EF\u25F8-\u25FC\u2605\u2606\u260E\u2640\u2642\u2660\u2663\u2665\u2666\u266A\u266D-\u266F\u2713\u2717\u2720\u2736\u2758\u2772\u2773\u27C8\u27C9\u27E6-\u27ED\u27F5-\u27FA\u27FC\u27FF\u2902-\u2905\u290C-\u2913\u2916\u2919-\u2920\u2923-\u292A\u2933\u2935-\u2939\u293C\u293D\u2945\u2948-\u294B\u294E-\u2976\u2978\u2979\u297B-\u297F\u2985\u2986\u298B-\u2996\u299A\u299C\u299D\u29A4-\u29B7\u29B9\u29BB\u29BC\u29BE-\u29C5\u29C9\u29CD-\u29D0\u29DC-\u29DE\u29E3-\u29E5\u29EB\u29F4\u29F6\u2A00-\u2A02\u2A04\u2A06\u2A0C\u2A0D\u2A10-\u2A17\u2A22-\u2A27\u2A29\u2A2A\u2A2D-\u2A31\u2A33-\u2A3C\u2A3F\u2A40\u2A42-\u2A4D\u2A50\u2A53-\u2A58\u2A5A-\u2A5D\u2A5F\u2A66\u2A6A\u2A6D-\u2A75\u2A77-\u2A9A\u2A9D-\u2AA2\u2AA4-\u2AB0\u2AB3-\u2AC8\u2ACB\u2ACC\u2ACF-\u2ADB\u2AE4\u2AE6-\u2AE9\u2AEB-\u2AF3\u2AFD\uFB00-\uFB04]|\uD835[\uDC9C\uDC9E\uDC9F\uDCA2\uDCA5\uDCA6\uDCA9-\uDCAC\uDCAE-\uDCB9\uDCBB\uDCBD-\uDCC3\uDCC5-\uDCCF\uDD04\uDD05\uDD07-\uDD0A\uDD0D-\uDD14\uDD16-\uDD1C\uDD1E-\uDD39\uDD3B-\uDD3E\uDD40-\uDD44\uDD46\uDD4A-\uDD50\uDD52-\uDD6B]/g;
		var encodeMap = {'\xAD':'shy','\u200C':'zwnj','\u200D':'zwj','\u200E':'lrm','\u2063':'ic','\u2062':'it','\u2061':'af','\u200F':'rlm','\u200B':'ZeroWidthSpace','\u2060':'NoBreak','\u0311':'DownBreve','\u20DB':'tdot','\u20DC':'DotDot','\t':'Tab','\n':'NewLine','\u2008':'puncsp','\u205F':'MediumSpace','\u2009':'thinsp','\u200A':'hairsp','\u2004':'emsp13','\u2002':'ensp','\u2005':'emsp14','\u2003':'emsp','\u2007':'numsp','\xA0':'nbsp','\u205F\u200A':'ThickSpace','\u203E':'oline','_':'lowbar','\u2010':'dash','\u2013':'ndash','\u2014':'mdash','\u2015':'horbar',',':'comma',';':'semi','\u204F':'bsemi',':':'colon','\u2A74':'Colone','!':'excl','\xA1':'iexcl','?':'quest','\xBF':'iquest','.':'period','\u2025':'nldr','\u2026':'mldr','\xB7':'middot','\'':'apos','\u2018':'lsquo','\u2019':'rsquo','\u201A':'sbquo','\u2039':'lsaquo','\u203A':'rsaquo','"':'quot','\u201C':'ldquo','\u201D':'rdquo','\u201E':'bdquo','\xAB':'laquo','\xBB':'raquo','(':'lpar',')':'rpar','[':'lsqb',']':'rsqb','{':'lcub','}':'rcub','\u2308':'lceil','\u2309':'rceil','\u230A':'lfloor','\u230B':'rfloor','\u2985':'lopar','\u2986':'ropar','\u298B':'lbrke','\u298C':'rbrke','\u298D':'lbrkslu','\u298E':'rbrksld','\u298F':'lbrksld','\u2990':'rbrkslu','\u2991':'langd','\u2992':'rangd','\u2993':'lparlt','\u2994':'rpargt','\u2995':'gtlPar','\u2996':'ltrPar','\u27E6':'lobrk','\u27E7':'robrk','\u27E8':'lang','\u27E9':'rang','\u27EA':'Lang','\u27EB':'Rang','\u27EC':'loang','\u27ED':'roang','\u2772':'lbbrk','\u2773':'rbbrk','\u2016':'Vert','\xA7':'sect','\xB6':'para','@':'commat','*':'ast','/':'sol','undefined':null,'&':'amp','#':'num','%':'percnt','\u2030':'permil','\u2031':'pertenk','\u2020':'dagger','\u2021':'Dagger','\u2022':'bull','\u2043':'hybull','\u2032':'prime','\u2033':'Prime','\u2034':'tprime','\u2057':'qprime','\u2035':'bprime','\u2041':'caret','`':'grave','\xB4':'acute','\u02DC':'tilde','^':'Hat','\xAF':'macr','\u02D8':'breve','\u02D9':'dot','\xA8':'die','\u02DA':'ring','\u02DD':'dblac','\xB8':'cedil','\u02DB':'ogon','\u02C6':'circ','\u02C7':'caron','\xB0':'deg','\xA9':'copy','\xAE':'reg','\u2117':'copysr','\u2118':'wp','\u211E':'rx','\u2127':'mho','\u2129':'iiota','\u2190':'larr','\u219A':'nlarr','\u2192':'rarr','\u219B':'nrarr','\u2191':'uarr','\u2193':'darr','\u2194':'harr','\u21AE':'nharr','\u2195':'varr','\u2196':'nwarr','\u2197':'nearr','\u2198':'searr','\u2199':'swarr','\u219D':'rarrw','\u219D\u0338':'nrarrw','\u219E':'Larr','\u219F':'Uarr','\u21A0':'Rarr','\u21A1':'Darr','\u21A2':'larrtl','\u21A3':'rarrtl','\u21A4':'mapstoleft','\u21A5':'mapstoup','\u21A6':'map','\u21A7':'mapstodown','\u21A9':'larrhk','\u21AA':'rarrhk','\u21AB':'larrlp','\u21AC':'rarrlp','\u21AD':'harrw','\u21B0':'lsh','\u21B1':'rsh','\u21B2':'ldsh','\u21B3':'rdsh','\u21B5':'crarr','\u21B6':'cularr','\u21B7':'curarr','\u21BA':'olarr','\u21BB':'orarr','\u21BC':'lharu','\u21BD':'lhard','\u21BE':'uharr','\u21BF':'uharl','\u21C0':'rharu','\u21C1':'rhard','\u21C2':'dharr','\u21C3':'dharl','\u21C4':'rlarr','\u21C5':'udarr','\u21C6':'lrarr','\u21C7':'llarr','\u21C8':'uuarr','\u21C9':'rrarr','\u21CA':'ddarr','\u21CB':'lrhar','\u21CC':'rlhar','\u21D0':'lArr','\u21CD':'nlArr','\u21D1':'uArr','\u21D2':'rArr','\u21CF':'nrArr','\u21D3':'dArr','\u21D4':'iff','\u21CE':'nhArr','\u21D5':'vArr','\u21D6':'nwArr','\u21D7':'neArr','\u21D8':'seArr','\u21D9':'swArr','\u21DA':'lAarr','\u21DB':'rAarr','\u21DD':'zigrarr','\u21E4':'larrb','\u21E5':'rarrb','\u21F5':'duarr','\u21FD':'loarr','\u21FE':'roarr','\u21FF':'hoarr','\u2200':'forall','\u2201':'comp','\u2202':'part','\u2202\u0338':'npart','\u2203':'exist','\u2204':'nexist','\u2205':'empty','\u2207':'Del','\u2208':'in','\u2209':'notin','\u220B':'ni','\u220C':'notni','\u03F6':'bepsi','\u220F':'prod','\u2210':'coprod','\u2211':'sum','+':'plus','\xB1':'pm','\xF7':'div','\xD7':'times','<':'lt','\u226E':'nlt','<\u20D2':'nvlt','=':'equals','\u2260':'ne','=\u20E5':'bne','\u2A75':'Equal','>':'gt','\u226F':'ngt','>\u20D2':'nvgt','\xAC':'not','|':'vert','\xA6':'brvbar','\u2212':'minus','\u2213':'mp','\u2214':'plusdo','\u2044':'frasl','\u2216':'setmn','\u2217':'lowast','\u2218':'compfn','\u221A':'Sqrt','\u221D':'prop','\u221E':'infin','\u221F':'angrt','\u2220':'ang','\u2220\u20D2':'nang','\u2221':'angmsd','\u2222':'angsph','\u2223':'mid','\u2224':'nmid','\u2225':'par','\u2226':'npar','\u2227':'and','\u2228':'or','\u2229':'cap','\u2229\uFE00':'caps','\u222A':'cup','\u222A\uFE00':'cups','\u222B':'int','\u222C':'Int','\u222D':'tint','\u2A0C':'qint','\u222E':'oint','\u222F':'Conint','\u2230':'Cconint','\u2231':'cwint','\u2232':'cwconint','\u2233':'awconint','\u2234':'there4','\u2235':'becaus','\u2236':'ratio','\u2237':'Colon','\u2238':'minusd','\u223A':'mDDot','\u223B':'homtht','\u223C':'sim','\u2241':'nsim','\u223C\u20D2':'nvsim','\u223D':'bsim','\u223D\u0331':'race','\u223E':'ac','\u223E\u0333':'acE','\u223F':'acd','\u2240':'wr','\u2242':'esim','\u2242\u0338':'nesim','\u2243':'sime','\u2244':'nsime','\u2245':'cong','\u2247':'ncong','\u2246':'simne','\u2248':'ap','\u2249':'nap','\u224A':'ape','\u224B':'apid','\u224B\u0338':'napid','\u224C':'bcong','\u224D':'CupCap','\u226D':'NotCupCap','\u224D\u20D2':'nvap','\u224E':'bump','\u224E\u0338':'nbump','\u224F':'bumpe','\u224F\u0338':'nbumpe','\u2250':'doteq','\u2250\u0338':'nedot','\u2251':'eDot','\u2252':'efDot','\u2253':'erDot','\u2254':'colone','\u2255':'ecolon','\u2256':'ecir','\u2257':'cire','\u2259':'wedgeq','\u225A':'veeeq','\u225C':'trie','\u225F':'equest','\u2261':'equiv','\u2262':'nequiv','\u2261\u20E5':'bnequiv','\u2264':'le','\u2270':'nle','\u2264\u20D2':'nvle','\u2265':'ge','\u2271':'nge','\u2265\u20D2':'nvge','\u2266':'lE','\u2266\u0338':'nlE','\u2267':'gE','\u2267\u0338':'ngE','\u2268\uFE00':'lvnE','\u2268':'lnE','\u2269':'gnE','\u2269\uFE00':'gvnE','\u226A':'ll','\u226A\u0338':'nLtv','\u226A\u20D2':'nLt','\u226B':'gg','\u226B\u0338':'nGtv','\u226B\u20D2':'nGt','\u226C':'twixt','\u2272':'lsim','\u2274':'nlsim','\u2273':'gsim','\u2275':'ngsim','\u2276':'lg','\u2278':'ntlg','\u2277':'gl','\u2279':'ntgl','\u227A':'pr','\u2280':'npr','\u227B':'sc','\u2281':'nsc','\u227C':'prcue','\u22E0':'nprcue','\u227D':'sccue','\u22E1':'nsccue','\u227E':'prsim','\u227F':'scsim','\u227F\u0338':'NotSucceedsTilde','\u2282':'sub','\u2284':'nsub','\u2282\u20D2':'vnsub','\u2283':'sup','\u2285':'nsup','\u2283\u20D2':'vnsup','\u2286':'sube','\u2288':'nsube','\u2287':'supe','\u2289':'nsupe','\u228A\uFE00':'vsubne','\u228A':'subne','\u228B\uFE00':'vsupne','\u228B':'supne','\u228D':'cupdot','\u228E':'uplus','\u228F':'sqsub','\u228F\u0338':'NotSquareSubset','\u2290':'sqsup','\u2290\u0338':'NotSquareSuperset','\u2291':'sqsube','\u22E2':'nsqsube','\u2292':'sqsupe','\u22E3':'nsqsupe','\u2293':'sqcap','\u2293\uFE00':'sqcaps','\u2294':'sqcup','\u2294\uFE00':'sqcups','\u2295':'oplus','\u2296':'ominus','\u2297':'otimes','\u2298':'osol','\u2299':'odot','\u229A':'ocir','\u229B':'oast','\u229D':'odash','\u229E':'plusb','\u229F':'minusb','\u22A0':'timesb','\u22A1':'sdotb','\u22A2':'vdash','\u22AC':'nvdash','\u22A3':'dashv','\u22A4':'top','\u22A5':'bot','\u22A7':'models','\u22A8':'vDash','\u22AD':'nvDash','\u22A9':'Vdash','\u22AE':'nVdash','\u22AA':'Vvdash','\u22AB':'VDash','\u22AF':'nVDash','\u22B0':'prurel','\u22B2':'vltri','\u22EA':'nltri','\u22B3':'vrtri','\u22EB':'nrtri','\u22B4':'ltrie','\u22EC':'nltrie','\u22B4\u20D2':'nvltrie','\u22B5':'rtrie','\u22ED':'nrtrie','\u22B5\u20D2':'nvrtrie','\u22B6':'origof','\u22B7':'imof','\u22B8':'mumap','\u22B9':'hercon','\u22BA':'intcal','\u22BB':'veebar','\u22BD':'barvee','\u22BE':'angrtvb','\u22BF':'lrtri','\u22C0':'Wedge','\u22C1':'Vee','\u22C2':'xcap','\u22C3':'xcup','\u22C4':'diam','\u22C5':'sdot','\u22C6':'Star','\u22C7':'divonx','\u22C8':'bowtie','\u22C9':'ltimes','\u22CA':'rtimes','\u22CB':'lthree','\u22CC':'rthree','\u22CD':'bsime','\u22CE':'cuvee','\u22CF':'cuwed','\u22D0':'Sub','\u22D1':'Sup','\u22D2':'Cap','\u22D3':'Cup','\u22D4':'fork','\u22D5':'epar','\u22D6':'ltdot','\u22D7':'gtdot','\u22D8':'Ll','\u22D8\u0338':'nLl','\u22D9':'Gg','\u22D9\u0338':'nGg','\u22DA\uFE00':'lesg','\u22DA':'leg','\u22DB':'gel','\u22DB\uFE00':'gesl','\u22DE':'cuepr','\u22DF':'cuesc','\u22E6':'lnsim','\u22E7':'gnsim','\u22E8':'prnsim','\u22E9':'scnsim','\u22EE':'vellip','\u22EF':'ctdot','\u22F0':'utdot','\u22F1':'dtdot','\u22F2':'disin','\u22F3':'isinsv','\u22F4':'isins','\u22F5':'isindot','\u22F5\u0338':'notindot','\u22F6':'notinvc','\u22F7':'notinvb','\u22F9':'isinE','\u22F9\u0338':'notinE','\u22FA':'nisd','\u22FB':'xnis','\u22FC':'nis','\u22FD':'notnivc','\u22FE':'notnivb','\u2305':'barwed','\u2306':'Barwed','\u230C':'drcrop','\u230D':'dlcrop','\u230E':'urcrop','\u230F':'ulcrop','\u2310':'bnot','\u2312':'profline','\u2313':'profsurf','\u2315':'telrec','\u2316':'target','\u231C':'ulcorn','\u231D':'urcorn','\u231E':'dlcorn','\u231F':'drcorn','\u2322':'frown','\u2323':'smile','\u232D':'cylcty','\u232E':'profalar','\u2336':'topbot','\u233D':'ovbar','\u233F':'solbar','\u237C':'angzarr','\u23B0':'lmoust','\u23B1':'rmoust','\u23B4':'tbrk','\u23B5':'bbrk','\u23B6':'bbrktbrk','\u23DC':'OverParenthesis','\u23DD':'UnderParenthesis','\u23DE':'OverBrace','\u23DF':'UnderBrace','\u23E2':'trpezium','\u23E7':'elinters','\u2423':'blank','\u2500':'boxh','\u2502':'boxv','\u250C':'boxdr','\u2510':'boxdl','\u2514':'boxur','\u2518':'boxul','\u251C':'boxvr','\u2524':'boxvl','\u252C':'boxhd','\u2534':'boxhu','\u253C':'boxvh','\u2550':'boxH','\u2551':'boxV','\u2552':'boxdR','\u2553':'boxDr','\u2554':'boxDR','\u2555':'boxdL','\u2556':'boxDl','\u2557':'boxDL','\u2558':'boxuR','\u2559':'boxUr','\u255A':'boxUR','\u255B':'boxuL','\u255C':'boxUl','\u255D':'boxUL','\u255E':'boxvR','\u255F':'boxVr','\u2560':'boxVR','\u2561':'boxvL','\u2562':'boxVl','\u2563':'boxVL','\u2564':'boxHd','\u2565':'boxhD','\u2566':'boxHD','\u2567':'boxHu','\u2568':'boxhU','\u2569':'boxHU','\u256A':'boxvH','\u256B':'boxVh','\u256C':'boxVH','\u2580':'uhblk','\u2584':'lhblk','\u2588':'block','\u2591':'blk14','\u2592':'blk12','\u2593':'blk34','\u25A1':'squ','\u25AA':'squf','\u25AB':'EmptyVerySmallSquare','\u25AD':'rect','\u25AE':'marker','\u25B1':'fltns','\u25B3':'xutri','\u25B4':'utrif','\u25B5':'utri','\u25B8':'rtrif','\u25B9':'rtri','\u25BD':'xdtri','\u25BE':'dtrif','\u25BF':'dtri','\u25C2':'ltrif','\u25C3':'ltri','\u25CA':'loz','\u25CB':'cir','\u25EC':'tridot','\u25EF':'xcirc','\u25F8':'ultri','\u25F9':'urtri','\u25FA':'lltri','\u25FB':'EmptySmallSquare','\u25FC':'FilledSmallSquare','\u2605':'starf','\u2606':'star','\u260E':'phone','\u2640':'female','\u2642':'male','\u2660':'spades','\u2663':'clubs','\u2665':'hearts','\u2666':'diams','\u266A':'sung','\u2713':'check','\u2717':'cross','\u2720':'malt','\u2736':'sext','\u2758':'VerticalSeparator','\u27C8':'bsolhsub','\u27C9':'suphsol','\u27F5':'xlarr','\u27F6':'xrarr','\u27F7':'xharr','\u27F8':'xlArr','\u27F9':'xrArr','\u27FA':'xhArr','\u27FC':'xmap','\u27FF':'dzigrarr','\u2902':'nvlArr','\u2903':'nvrArr','\u2904':'nvHarr','\u2905':'Map','\u290C':'lbarr','\u290D':'rbarr','\u290E':'lBarr','\u290F':'rBarr','\u2910':'RBarr','\u2911':'DDotrahd','\u2912':'UpArrowBar','\u2913':'DownArrowBar','\u2916':'Rarrtl','\u2919':'latail','\u291A':'ratail','\u291B':'lAtail','\u291C':'rAtail','\u291D':'larrfs','\u291E':'rarrfs','\u291F':'larrbfs','\u2920':'rarrbfs','\u2923':'nwarhk','\u2924':'nearhk','\u2925':'searhk','\u2926':'swarhk','\u2927':'nwnear','\u2928':'toea','\u2929':'tosa','\u292A':'swnwar','\u2933':'rarrc','\u2933\u0338':'nrarrc','\u2935':'cudarrr','\u2936':'ldca','\u2937':'rdca','\u2938':'cudarrl','\u2939':'larrpl','\u293C':'curarrm','\u293D':'cularrp','\u2945':'rarrpl','\u2948':'harrcir','\u2949':'Uarrocir','\u294A':'lurdshar','\u294B':'ldrushar','\u294E':'LeftRightVector','\u294F':'RightUpDownVector','\u2950':'DownLeftRightVector','\u2951':'LeftUpDownVector','\u2952':'LeftVectorBar','\u2953':'RightVectorBar','\u2954':'RightUpVectorBar','\u2955':'RightDownVectorBar','\u2956':'DownLeftVectorBar','\u2957':'DownRightVectorBar','\u2958':'LeftUpVectorBar','\u2959':'LeftDownVectorBar','\u295A':'LeftTeeVector','\u295B':'RightTeeVector','\u295C':'RightUpTeeVector','\u295D':'RightDownTeeVector','\u295E':'DownLeftTeeVector','\u295F':'DownRightTeeVector','\u2960':'LeftUpTeeVector','\u2961':'LeftDownTeeVector','\u2962':'lHar','\u2963':'uHar','\u2964':'rHar','\u2965':'dHar','\u2966':'luruhar','\u2967':'ldrdhar','\u2968':'ruluhar','\u2969':'rdldhar','\u296A':'lharul','\u296B':'llhard','\u296C':'rharul','\u296D':'lrhard','\u296E':'udhar','\u296F':'duhar','\u2970':'RoundImplies','\u2971':'erarr','\u2972':'simrarr','\u2973':'larrsim','\u2974':'rarrsim','\u2975':'rarrap','\u2976':'ltlarr','\u2978':'gtrarr','\u2979':'subrarr','\u297B':'suplarr','\u297C':'lfisht','\u297D':'rfisht','\u297E':'ufisht','\u297F':'dfisht','\u299A':'vzigzag','\u299C':'vangrt','\u299D':'angrtvbd','\u29A4':'ange','\u29A5':'range','\u29A6':'dwangle','\u29A7':'uwangle','\u29A8':'angmsdaa','\u29A9':'angmsdab','\u29AA':'angmsdac','\u29AB':'angmsdad','\u29AC':'angmsdae','\u29AD':'angmsdaf','\u29AE':'angmsdag','\u29AF':'angmsdah','\u29B0':'bemptyv','\u29B1':'demptyv','\u29B2':'cemptyv','\u29B3':'raemptyv','\u29B4':'laemptyv','\u29B5':'ohbar','\u29B6':'omid','\u29B7':'opar','\u29B9':'operp','\u29BB':'olcross','\u29BC':'odsold','\u29BE':'olcir','\u29BF':'ofcir','\u29C0':'olt','\u29C1':'ogt','\u29C2':'cirscir','\u29C3':'cirE','\u29C4':'solb','\u29C5':'bsolb','\u29C9':'boxbox','\u29CD':'trisb','\u29CE':'rtriltri','\u29CF':'LeftTriangleBar','\u29CF\u0338':'NotLeftTriangleBar','\u29D0':'RightTriangleBar','\u29D0\u0338':'NotRightTriangleBar','\u29DC':'iinfin','\u29DD':'infintie','\u29DE':'nvinfin','\u29E3':'eparsl','\u29E4':'smeparsl','\u29E5':'eqvparsl','\u29EB':'lozf','\u29F4':'RuleDelayed','\u29F6':'dsol','\u2A00':'xodot','\u2A01':'xoplus','\u2A02':'xotime','\u2A04':'xuplus','\u2A06':'xsqcup','\u2A0D':'fpartint','\u2A10':'cirfnint','\u2A11':'awint','\u2A12':'rppolint','\u2A13':'scpolint','\u2A14':'npolint','\u2A15':'pointint','\u2A16':'quatint','\u2A17':'intlarhk','\u2A22':'pluscir','\u2A23':'plusacir','\u2A24':'simplus','\u2A25':'plusdu','\u2A26':'plussim','\u2A27':'plustwo','\u2A29':'mcomma','\u2A2A':'minusdu','\u2A2D':'loplus','\u2A2E':'roplus','\u2A2F':'Cross','\u2A30':'timesd','\u2A31':'timesbar','\u2A33':'smashp','\u2A34':'lotimes','\u2A35':'rotimes','\u2A36':'otimesas','\u2A37':'Otimes','\u2A38':'odiv','\u2A39':'triplus','\u2A3A':'triminus','\u2A3B':'tritime','\u2A3C':'iprod','\u2A3F':'amalg','\u2A40':'capdot','\u2A42':'ncup','\u2A43':'ncap','\u2A44':'capand','\u2A45':'cupor','\u2A46':'cupcap','\u2A47':'capcup','\u2A48':'cupbrcap','\u2A49':'capbrcup','\u2A4A':'cupcup','\u2A4B':'capcap','\u2A4C':'ccups','\u2A4D':'ccaps','\u2A50':'ccupssm','\u2A53':'And','\u2A54':'Or','\u2A55':'andand','\u2A56':'oror','\u2A57':'orslope','\u2A58':'andslope','\u2A5A':'andv','\u2A5B':'orv','\u2A5C':'andd','\u2A5D':'ord','\u2A5F':'wedbar','\u2A66':'sdote','\u2A6A':'simdot','\u2A6D':'congdot','\u2A6D\u0338':'ncongdot','\u2A6E':'easter','\u2A6F':'apacir','\u2A70':'apE','\u2A70\u0338':'napE','\u2A71':'eplus','\u2A72':'pluse','\u2A73':'Esim','\u2A77':'eDDot','\u2A78':'equivDD','\u2A79':'ltcir','\u2A7A':'gtcir','\u2A7B':'ltquest','\u2A7C':'gtquest','\u2A7D':'les','\u2A7D\u0338':'nles','\u2A7E':'ges','\u2A7E\u0338':'nges','\u2A7F':'lesdot','\u2A80':'gesdot','\u2A81':'lesdoto','\u2A82':'gesdoto','\u2A83':'lesdotor','\u2A84':'gesdotol','\u2A85':'lap','\u2A86':'gap','\u2A87':'lne','\u2A88':'gne','\u2A89':'lnap','\u2A8A':'gnap','\u2A8B':'lEg','\u2A8C':'gEl','\u2A8D':'lsime','\u2A8E':'gsime','\u2A8F':'lsimg','\u2A90':'gsiml','\u2A91':'lgE','\u2A92':'glE','\u2A93':'lesges','\u2A94':'gesles','\u2A95':'els','\u2A96':'egs','\u2A97':'elsdot','\u2A98':'egsdot','\u2A99':'el','\u2A9A':'eg','\u2A9D':'siml','\u2A9E':'simg','\u2A9F':'simlE','\u2AA0':'simgE','\u2AA1':'LessLess','\u2AA1\u0338':'NotNestedLessLess','\u2AA2':'GreaterGreater','\u2AA2\u0338':'NotNestedGreaterGreater','\u2AA4':'glj','\u2AA5':'gla','\u2AA6':'ltcc','\u2AA7':'gtcc','\u2AA8':'lescc','\u2AA9':'gescc','\u2AAA':'smt','\u2AAB':'lat','\u2AAC':'smte','\u2AAC\uFE00':'smtes','\u2AAD':'late','\u2AAD\uFE00':'lates','\u2AAE':'bumpE','\u2AAF':'pre','\u2AAF\u0338':'npre','\u2AB0':'sce','\u2AB0\u0338':'nsce','\u2AB3':'prE','\u2AB4':'scE','\u2AB5':'prnE','\u2AB6':'scnE','\u2AB7':'prap','\u2AB8':'scap','\u2AB9':'prnap','\u2ABA':'scnap','\u2ABB':'Pr','\u2ABC':'Sc','\u2ABD':'subdot','\u2ABE':'supdot','\u2ABF':'subplus','\u2AC0':'supplus','\u2AC1':'submult','\u2AC2':'supmult','\u2AC3':'subedot','\u2AC4':'supedot','\u2AC5':'subE','\u2AC5\u0338':'nsubE','\u2AC6':'supE','\u2AC6\u0338':'nsupE','\u2AC7':'subsim','\u2AC8':'supsim','\u2ACB\uFE00':'vsubnE','\u2ACB':'subnE','\u2ACC\uFE00':'vsupnE','\u2ACC':'supnE','\u2ACF':'csub','\u2AD0':'csup','\u2AD1':'csube','\u2AD2':'csupe','\u2AD3':'subsup','\u2AD4':'supsub','\u2AD5':'subsub','\u2AD6':'supsup','\u2AD7':'suphsub','\u2AD8':'supdsub','\u2AD9':'forkv','\u2ADA':'topfork','\u2ADB':'mlcp','\u2AE4':'Dashv','\u2AE6':'Vdashl','\u2AE7':'Barv','\u2AE8':'vBar','\u2AE9':'vBarv','\u2AEB':'Vbar','\u2AEC':'Not','\u2AED':'bNot','\u2AEE':'rnmid','\u2AEF':'cirmid','\u2AF0':'midcir','\u2AF1':'topcir','\u2AF2':'nhpar','\u2AF3':'parsim','\u2AFD':'parsl','\u2AFD\u20E5':'nparsl','\u266D':'flat','\u266E':'natur','\u266F':'sharp','\xA4':'curren','\xA2':'cent','$':'dollar','\xA3':'pound','\xA5':'yen','\u20AC':'euro','\xB9':'sup1','\xBD':'half','\u2153':'frac13','\xBC':'frac14','\u2155':'frac15','\u2159':'frac16','\u215B':'frac18','\xB2':'sup2','\u2154':'frac23','\u2156':'frac25','\xB3':'sup3','\xBE':'frac34','\u2157':'frac35','\u215C':'frac38','\u2158':'frac45','\u215A':'frac56','\u215D':'frac58','\u215E':'frac78','\uD835\uDCB6':'ascr','\uD835\uDD52':'aopf','\uD835\uDD1E':'afr','\uD835\uDD38':'Aopf','\uD835\uDD04':'Afr','\uD835\uDC9C':'Ascr','\xAA':'ordf','\xE1':'aacute','\xC1':'Aacute','\xE0':'agrave','\xC0':'Agrave','\u0103':'abreve','\u0102':'Abreve','\xE2':'acirc','\xC2':'Acirc','\xE5':'aring','\xC5':'angst','\xE4':'auml','\xC4':'Auml','\xE3':'atilde','\xC3':'Atilde','\u0105':'aogon','\u0104':'Aogon','\u0101':'amacr','\u0100':'Amacr','\xE6':'aelig','\xC6':'AElig','\uD835\uDCB7':'bscr','\uD835\uDD53':'bopf','\uD835\uDD1F':'bfr','\uD835\uDD39':'Bopf','\u212C':'Bscr','\uD835\uDD05':'Bfr','\uD835\uDD20':'cfr','\uD835\uDCB8':'cscr','\uD835\uDD54':'copf','\u212D':'Cfr','\uD835\uDC9E':'Cscr','\u2102':'Copf','\u0107':'cacute','\u0106':'Cacute','\u0109':'ccirc','\u0108':'Ccirc','\u010D':'ccaron','\u010C':'Ccaron','\u010B':'cdot','\u010A':'Cdot','\xE7':'ccedil','\xC7':'Ccedil','\u2105':'incare','\uD835\uDD21':'dfr','\u2146':'dd','\uD835\uDD55':'dopf','\uD835\uDCB9':'dscr','\uD835\uDC9F':'Dscr','\uD835\uDD07':'Dfr','\u2145':'DD','\uD835\uDD3B':'Dopf','\u010F':'dcaron','\u010E':'Dcaron','\u0111':'dstrok','\u0110':'Dstrok','\xF0':'eth','\xD0':'ETH','\u2147':'ee','\u212F':'escr','\uD835\uDD22':'efr','\uD835\uDD56':'eopf','\u2130':'Escr','\uD835\uDD08':'Efr','\uD835\uDD3C':'Eopf','\xE9':'eacute','\xC9':'Eacute','\xE8':'egrave','\xC8':'Egrave','\xEA':'ecirc','\xCA':'Ecirc','\u011B':'ecaron','\u011A':'Ecaron','\xEB':'euml','\xCB':'Euml','\u0117':'edot','\u0116':'Edot','\u0119':'eogon','\u0118':'Eogon','\u0113':'emacr','\u0112':'Emacr','\uD835\uDD23':'ffr','\uD835\uDD57':'fopf','\uD835\uDCBB':'fscr','\uD835\uDD09':'Ffr','\uD835\uDD3D':'Fopf','\u2131':'Fscr','\uFB00':'fflig','\uFB03':'ffilig','\uFB04':'ffllig','\uFB01':'filig','fj':'fjlig','\uFB02':'fllig','\u0192':'fnof','\u210A':'gscr','\uD835\uDD58':'gopf','\uD835\uDD24':'gfr','\uD835\uDCA2':'Gscr','\uD835\uDD3E':'Gopf','\uD835\uDD0A':'Gfr','\u01F5':'gacute','\u011F':'gbreve','\u011E':'Gbreve','\u011D':'gcirc','\u011C':'Gcirc','\u0121':'gdot','\u0120':'Gdot','\u0122':'Gcedil','\uD835\uDD25':'hfr','\u210E':'planckh','\uD835\uDCBD':'hscr','\uD835\uDD59':'hopf','\u210B':'Hscr','\u210C':'Hfr','\u210D':'Hopf','\u0125':'hcirc','\u0124':'Hcirc','\u210F':'hbar','\u0127':'hstrok','\u0126':'Hstrok','\uD835\uDD5A':'iopf','\uD835\uDD26':'ifr','\uD835\uDCBE':'iscr','\u2148':'ii','\uD835\uDD40':'Iopf','\u2110':'Iscr','\u2111':'Im','\xED':'iacute','\xCD':'Iacute','\xEC':'igrave','\xCC':'Igrave','\xEE':'icirc','\xCE':'Icirc','\xEF':'iuml','\xCF':'Iuml','\u0129':'itilde','\u0128':'Itilde','\u0130':'Idot','\u012F':'iogon','\u012E':'Iogon','\u012B':'imacr','\u012A':'Imacr','\u0133':'ijlig','\u0132':'IJlig','\u0131':'imath','\uD835\uDCBF':'jscr','\uD835\uDD5B':'jopf','\uD835\uDD27':'jfr','\uD835\uDCA5':'Jscr','\uD835\uDD0D':'Jfr','\uD835\uDD41':'Jopf','\u0135':'jcirc','\u0134':'Jcirc','\u0237':'jmath','\uD835\uDD5C':'kopf','\uD835\uDCC0':'kscr','\uD835\uDD28':'kfr','\uD835\uDCA6':'Kscr','\uD835\uDD42':'Kopf','\uD835\uDD0E':'Kfr','\u0137':'kcedil','\u0136':'Kcedil','\uD835\uDD29':'lfr','\uD835\uDCC1':'lscr','\u2113':'ell','\uD835\uDD5D':'lopf','\u2112':'Lscr','\uD835\uDD0F':'Lfr','\uD835\uDD43':'Lopf','\u013A':'lacute','\u0139':'Lacute','\u013E':'lcaron','\u013D':'Lcaron','\u013C':'lcedil','\u013B':'Lcedil','\u0142':'lstrok','\u0141':'Lstrok','\u0140':'lmidot','\u013F':'Lmidot','\uD835\uDD2A':'mfr','\uD835\uDD5E':'mopf','\uD835\uDCC2':'mscr','\uD835\uDD10':'Mfr','\uD835\uDD44':'Mopf','\u2133':'Mscr','\uD835\uDD2B':'nfr','\uD835\uDD5F':'nopf','\uD835\uDCC3':'nscr','\u2115':'Nopf','\uD835\uDCA9':'Nscr','\uD835\uDD11':'Nfr','\u0144':'nacute','\u0143':'Nacute','\u0148':'ncaron','\u0147':'Ncaron','\xF1':'ntilde','\xD1':'Ntilde','\u0146':'ncedil','\u0145':'Ncedil','\u2116':'numero','\u014B':'eng','\u014A':'ENG','\uD835\uDD60':'oopf','\uD835\uDD2C':'ofr','\u2134':'oscr','\uD835\uDCAA':'Oscr','\uD835\uDD12':'Ofr','\uD835\uDD46':'Oopf','\xBA':'ordm','\xF3':'oacute','\xD3':'Oacute','\xF2':'ograve','\xD2':'Ograve','\xF4':'ocirc','\xD4':'Ocirc','\xF6':'ouml','\xD6':'Ouml','\u0151':'odblac','\u0150':'Odblac','\xF5':'otilde','\xD5':'Otilde','\xF8':'oslash','\xD8':'Oslash','\u014D':'omacr','\u014C':'Omacr','\u0153':'oelig','\u0152':'OElig','\uD835\uDD2D':'pfr','\uD835\uDCC5':'pscr','\uD835\uDD61':'popf','\u2119':'Popf','\uD835\uDD13':'Pfr','\uD835\uDCAB':'Pscr','\uD835\uDD62':'qopf','\uD835\uDD2E':'qfr','\uD835\uDCC6':'qscr','\uD835\uDCAC':'Qscr','\uD835\uDD14':'Qfr','\u211A':'Qopf','\u0138':'kgreen','\uD835\uDD2F':'rfr','\uD835\uDD63':'ropf','\uD835\uDCC7':'rscr','\u211B':'Rscr','\u211C':'Re','\u211D':'Ropf','\u0155':'racute','\u0154':'Racute','\u0159':'rcaron','\u0158':'Rcaron','\u0157':'rcedil','\u0156':'Rcedil','\uD835\uDD64':'sopf','\uD835\uDCC8':'sscr','\uD835\uDD30':'sfr','\uD835\uDD4A':'Sopf','\uD835\uDD16':'Sfr','\uD835\uDCAE':'Sscr','\u24C8':'oS','\u015B':'sacute','\u015A':'Sacute','\u015D':'scirc','\u015C':'Scirc','\u0161':'scaron','\u0160':'Scaron','\u015F':'scedil','\u015E':'Scedil','\xDF':'szlig','\uD835\uDD31':'tfr','\uD835\uDCC9':'tscr','\uD835\uDD65':'topf','\uD835\uDCAF':'Tscr','\uD835\uDD17':'Tfr','\uD835\uDD4B':'Topf','\u0165':'tcaron','\u0164':'Tcaron','\u0163':'tcedil','\u0162':'Tcedil','\u2122':'trade','\u0167':'tstrok','\u0166':'Tstrok','\uD835\uDCCA':'uscr','\uD835\uDD66':'uopf','\uD835\uDD32':'ufr','\uD835\uDD4C':'Uopf','\uD835\uDD18':'Ufr','\uD835\uDCB0':'Uscr','\xFA':'uacute','\xDA':'Uacute','\xF9':'ugrave','\xD9':'Ugrave','\u016D':'ubreve','\u016C':'Ubreve','\xFB':'ucirc','\xDB':'Ucirc','\u016F':'uring','\u016E':'Uring','\xFC':'uuml','\xDC':'Uuml','\u0171':'udblac','\u0170':'Udblac','\u0169':'utilde','\u0168':'Utilde','\u0173':'uogon','\u0172':'Uogon','\u016B':'umacr','\u016A':'Umacr','\uD835\uDD33':'vfr','\uD835\uDD67':'vopf','\uD835\uDCCB':'vscr','\uD835\uDD19':'Vfr','\uD835\uDD4D':'Vopf','\uD835\uDCB1':'Vscr','\uD835\uDD68':'wopf','\uD835\uDCCC':'wscr','\uD835\uDD34':'wfr','\uD835\uDCB2':'Wscr','\uD835\uDD4E':'Wopf','\uD835\uDD1A':'Wfr','\u0175':'wcirc','\u0174':'Wcirc','\uD835\uDD35':'xfr','\uD835\uDCCD':'xscr','\uD835\uDD69':'xopf','\uD835\uDD4F':'Xopf','\uD835\uDD1B':'Xfr','\uD835\uDCB3':'Xscr','\uD835\uDD36':'yfr','\uD835\uDCCE':'yscr','\uD835\uDD6A':'yopf','\uD835\uDCB4':'Yscr','\uD835\uDD1C':'Yfr','\uD835\uDD50':'Yopf','\xFD':'yacute','\xDD':'Yacute','\u0177':'ycirc','\u0176':'Ycirc','\xFF':'yuml','\u0178':'Yuml','\uD835\uDCCF':'zscr','\uD835\uDD37':'zfr','\uD835\uDD6B':'zopf','\u2128':'Zfr','\u2124':'Zopf','\uD835\uDCB5':'Zscr','\u017A':'zacute','\u0179':'Zacute','\u017E':'zcaron','\u017D':'Zcaron','\u017C':'zdot','\u017B':'Zdot','\u01B5':'imped','\xFE':'thorn','\xDE':'THORN','\u0149':'napos','\u03B1':'alpha','\u0391':'Alpha','\u03B2':'beta','\u0392':'Beta','\u03B3':'gamma','\u0393':'Gamma','\u03B4':'delta','\u0394':'Delta','\u03B5':'epsi','\u03F5':'epsiv','\u0395':'Epsilon','\u03DD':'gammad','\u03DC':'Gammad','\u03B6':'zeta','\u0396':'Zeta','\u03B7':'eta','\u0397':'Eta','\u03B8':'theta','\u03D1':'thetav','\u0398':'Theta','\u03B9':'iota','\u0399':'Iota','\u03BA':'kappa','\u03F0':'kappav','\u039A':'Kappa','\u03BB':'lambda','\u039B':'Lambda','\u03BC':'mu','\xB5':'micro','\u039C':'Mu','\u03BD':'nu','\u039D':'Nu','\u03BE':'xi','\u039E':'Xi','\u03BF':'omicron','\u039F':'Omicron','\u03C0':'pi','\u03D6':'piv','\u03A0':'Pi','\u03C1':'rho','\u03F1':'rhov','\u03A1':'Rho','\u03C3':'sigma','\u03A3':'Sigma','\u03C2':'sigmaf','\u03C4':'tau','\u03A4':'Tau','\u03C5':'upsi','\u03A5':'Upsilon','\u03D2':'Upsi','\u03C6':'phi','\u03D5':'phiv','\u03A6':'Phi','\u03C7':'chi','\u03A7':'Chi','\u03C8':'psi','\u03A8':'Psi','\u03C9':'omega','\u03A9':'ohm','\u0430':'acy','\u0410':'Acy','\u0431':'bcy','\u0411':'Bcy','\u0432':'vcy','\u0412':'Vcy','\u0433':'gcy','\u0413':'Gcy','\u0453':'gjcy','\u0403':'GJcy','\u0434':'dcy','\u0414':'Dcy','\u0452':'djcy','\u0402':'DJcy','\u0435':'iecy','\u0415':'IEcy','\u0451':'iocy','\u0401':'IOcy','\u0454':'jukcy','\u0404':'Jukcy','\u0436':'zhcy','\u0416':'ZHcy','\u0437':'zcy','\u0417':'Zcy','\u0455':'dscy','\u0405':'DScy','\u0438':'icy','\u0418':'Icy','\u0456':'iukcy','\u0406':'Iukcy','\u0457':'yicy','\u0407':'YIcy','\u0439':'jcy','\u0419':'Jcy','\u0458':'jsercy','\u0408':'Jsercy','\u043A':'kcy','\u041A':'Kcy','\u045C':'kjcy','\u040C':'KJcy','\u043B':'lcy','\u041B':'Lcy','\u0459':'ljcy','\u0409':'LJcy','\u043C':'mcy','\u041C':'Mcy','\u043D':'ncy','\u041D':'Ncy','\u045A':'njcy','\u040A':'NJcy','\u043E':'ocy','\u041E':'Ocy','\u043F':'pcy','\u041F':'Pcy','\u0440':'rcy','\u0420':'Rcy','\u0441':'scy','\u0421':'Scy','\u0442':'tcy','\u0422':'Tcy','\u045B':'tshcy','\u040B':'TSHcy','\u0443':'ucy','\u0423':'Ucy','\u045E':'ubrcy','\u040E':'Ubrcy','\u0444':'fcy','\u0424':'Fcy','\u0445':'khcy','\u0425':'KHcy','\u0446':'tscy','\u0426':'TScy','\u0447':'chcy','\u0427':'CHcy','\u045F':'dzcy','\u040F':'DZcy','\u0448':'shcy','\u0428':'SHcy','\u0449':'shchcy','\u0429':'SHCHcy','\u044A':'hardcy','\u042A':'HARDcy','\u044B':'ycy','\u042B':'Ycy','\u044C':'softcy','\u042C':'SOFTcy','\u044D':'ecy','\u042D':'Ecy','\u044E':'yucy','\u042E':'YUcy','\u044F':'yacy','\u042F':'YAcy','\u2135':'aleph','\u2136':'beth','\u2137':'gimel','\u2138':'daleth'};

		var regexEscape = /["&'<>`]/g;
		var escapeMap = {
			'"': '&quot;',
			'&': '&amp;',
			'\'': '&#x27;',
			'<': '&lt;',
			// See https://mathiasbynens.be/notes/ambiguous-ampersands: in HTML, the
			// following is not strictly necessary unless it’s part of a tag or an
			// unquoted attribute value. We’re only escaping it to support those
			// situations, and for XML support.
			'>': '&gt;',
			// In Internet Explorer ≤ 8, the backtick character can be used
			// to break out of (un)quoted attribute values or HTML comments.
			// See http://html5sec.org/#102, http://html5sec.org/#108, and
			// http://html5sec.org/#133.
			'`': '&#x60;'
		};

		var regexInvalidEntity = /&#(?:[xX][^a-fA-F0-9]|[^0-9xX])/;
		var regexInvalidRawCodePoint = /[\0-\x08\x0B\x0E-\x1F\x7F-\x9F\uFDD0-\uFDEF\uFFFE\uFFFF]|[\uD83F\uD87F\uD8BF\uD8FF\uD93F\uD97F\uD9BF\uD9FF\uDA3F\uDA7F\uDABF\uDAFF\uDB3F\uDB7F\uDBBF\uDBFF][\uDFFE\uDFFF]|[\uD800-\uDBFF](?![\uDC00-\uDFFF])|(?:[^\uD800-\uDBFF]|^)[\uDC00-\uDFFF]/;
		var regexDecode = /&(CounterClockwiseContourIntegral|DoubleLongLeftRightArrow|ClockwiseContourIntegral|NotNestedGreaterGreater|NotSquareSupersetEqual|DiacriticalDoubleAcute|NotRightTriangleEqual|NotSucceedsSlantEqual|NotPrecedesSlantEqual|CloseCurlyDoubleQuote|NegativeVeryThinSpace|DoubleContourIntegral|FilledVerySmallSquare|CapitalDifferentialD|OpenCurlyDoubleQuote|EmptyVerySmallSquare|NestedGreaterGreater|DoubleLongRightArrow|NotLeftTriangleEqual|NotGreaterSlantEqual|ReverseUpEquilibrium|DoubleLeftRightArrow|NotSquareSubsetEqual|NotDoubleVerticalBar|RightArrowLeftArrow|NotGreaterFullEqual|NotRightTriangleBar|SquareSupersetEqual|DownLeftRightVector|DoubleLongLeftArrow|leftrightsquigarrow|LeftArrowRightArrow|NegativeMediumSpace|blacktriangleright|RightDownVectorBar|PrecedesSlantEqual|RightDoubleBracket|SucceedsSlantEqual|NotLeftTriangleBar|RightTriangleEqual|SquareIntersection|RightDownTeeVector|ReverseEquilibrium|NegativeThickSpace|longleftrightarrow|Longleftrightarrow|LongLeftRightArrow|DownRightTeeVector|DownRightVectorBar|GreaterSlantEqual|SquareSubsetEqual|LeftDownVectorBar|LeftDoubleBracket|VerticalSeparator|rightleftharpoons|NotGreaterGreater|NotSquareSuperset|blacktriangleleft|blacktriangledown|NegativeThinSpace|LeftDownTeeVector|NotLessSlantEqual|leftrightharpoons|DoubleUpDownArrow|DoubleVerticalBar|LeftTriangleEqual|FilledSmallSquare|twoheadrightarrow|NotNestedLessLess|DownLeftTeeVector|DownLeftVectorBar|RightAngleBracket|NotTildeFullEqual|NotReverseElement|RightUpDownVector|DiacriticalTilde|NotSucceedsTilde|circlearrowright|NotPrecedesEqual|rightharpoondown|DoubleRightArrow|NotSucceedsEqual|NonBreakingSpace|NotRightTriangle|LessEqualGreater|RightUpTeeVector|LeftAngleBracket|GreaterFullEqual|DownArrowUpArrow|RightUpVectorBar|twoheadleftarrow|GreaterEqualLess|downharpoonright|RightTriangleBar|ntrianglerighteq|NotSupersetEqual|LeftUpDownVector|DiacriticalAcute|rightrightarrows|vartriangleright|UpArrowDownArrow|DiacriticalGrave|UnderParenthesis|EmptySmallSquare|LeftUpVectorBar|leftrightarrows|DownRightVector|downharpoonleft|trianglerighteq|ShortRightArrow|OverParenthesis|DoubleLeftArrow|DoubleDownArrow|NotSquareSubset|bigtriangledown|ntrianglelefteq|UpperRightArrow|curvearrowright|vartriangleleft|NotLeftTriangle|nleftrightarrow|LowerRightArrow|NotHumpDownHump|NotGreaterTilde|rightthreetimes|LeftUpTeeVector|NotGreaterEqual|straightepsilon|LeftTriangleBar|rightsquigarrow|ContourIntegral|rightleftarrows|CloseCurlyQuote|RightDownVector|LeftRightVector|nLeftrightarrow|leftharpoondown|circlearrowleft|SquareSuperset|OpenCurlyQuote|hookrightarrow|HorizontalLine|DiacriticalDot|NotLessGreater|ntriangleright|DoubleRightTee|InvisibleComma|InvisibleTimes|LowerLeftArrow|DownLeftVector|NotSubsetEqual|curvearrowleft|trianglelefteq|NotVerticalBar|TildeFullEqual|downdownarrows|NotGreaterLess|RightTeeVector|ZeroWidthSpace|looparrowright|LongRightArrow|doublebarwedge|ShortLeftArrow|ShortDownArrow|RightVectorBar|GreaterGreater|ReverseElement|rightharpoonup|LessSlantEqual|leftthreetimes|upharpoonright|rightarrowtail|LeftDownVector|Longrightarrow|NestedLessLess|UpperLeftArrow|nshortparallel|leftleftarrows|leftrightarrow|Leftrightarrow|LeftRightArrow|longrightarrow|upharpoonleft|RightArrowBar|ApplyFunction|LeftTeeVector|leftarrowtail|NotEqualTilde|varsubsetneqq|varsupsetneqq|RightTeeArrow|SucceedsEqual|SucceedsTilde|LeftVectorBar|SupersetEqual|hookleftarrow|DifferentialD|VerticalTilde|VeryThinSpace|blacktriangle|bigtriangleup|LessFullEqual|divideontimes|leftharpoonup|UpEquilibrium|ntriangleleft|RightTriangle|measuredangle|shortparallel|longleftarrow|Longleftarrow|LongLeftArrow|DoubleLeftTee|Poincareplane|PrecedesEqual|triangleright|DoubleUpArrow|RightUpVector|fallingdotseq|looparrowleft|PrecedesTilde|NotTildeEqual|NotTildeTilde|smallsetminus|Proportional|triangleleft|triangledown|UnderBracket|NotHumpEqual|exponentiale|ExponentialE|NotLessTilde|HilbertSpace|RightCeiling|blacklozenge|varsupsetneq|HumpDownHump|GreaterEqual|VerticalLine|LeftTeeArrow|NotLessEqual|DownTeeArrow|LeftTriangle|varsubsetneq|Intersection|NotCongruent|DownArrowBar|LeftUpVector|LeftArrowBar|risingdotseq|GreaterTilde|RoundImplies|SquareSubset|ShortUpArrow|NotSuperset|quaternions|precnapprox|backepsilon|preccurlyeq|OverBracket|blacksquare|MediumSpace|VerticalBar|circledcirc|circleddash|CircleMinus|CircleTimes|LessGreater|curlyeqprec|curlyeqsucc|diamondsuit|UpDownArrow|Updownarrow|RuleDelayed|Rrightarrow|updownarrow|RightVector|nRightarrow|nrightarrow|eqslantless|LeftCeiling|Equilibrium|SmallCircle|expectation|NotSucceeds|thickapprox|GreaterLess|SquareUnion|NotPrecedes|NotLessLess|straightphi|succnapprox|succcurlyeq|SubsetEqual|sqsupseteq|Proportion|Laplacetrf|ImaginaryI|supsetneqq|NotGreater|gtreqqless|NotElement|ThickSpace|TildeEqual|TildeTilde|Fouriertrf|rmoustache|EqualTilde|eqslantgtr|UnderBrace|LeftVector|UpArrowBar|nLeftarrow|nsubseteqq|subsetneqq|nsupseteqq|nleftarrow|succapprox|lessapprox|UpTeeArrow|upuparrows|curlywedge|lesseqqgtr|varepsilon|varnothing|RightFloor|complement|CirclePlus|sqsubseteq|Lleftarrow|circledast|RightArrow|Rightarrow|rightarrow|lmoustache|Bernoullis|precapprox|mapstoleft|mapstodown|longmapsto|dotsquare|downarrow|DoubleDot|nsubseteq|supsetneq|leftarrow|nsupseteq|subsetneq|ThinSpace|ngeqslant|subseteqq|HumpEqual|NotSubset|triangleq|NotCupCap|lesseqgtr|heartsuit|TripleDot|Leftarrow|Coproduct|Congruent|varpropto|complexes|gvertneqq|LeftArrow|LessTilde|supseteqq|MinusPlus|CircleDot|nleqslant|NotExists|gtreqless|nparallel|UnionPlus|LeftFloor|checkmark|CenterDot|centerdot|Mellintrf|gtrapprox|bigotimes|OverBrace|spadesuit|therefore|pitchfork|rationals|PlusMinus|Backslash|Therefore|DownBreve|backsimeq|backprime|DownArrow|nshortmid|Downarrow|lvertneqq|eqvparsl|imagline|imagpart|infintie|integers|Integral|intercal|LessLess|Uarrocir|intlarhk|sqsupset|angmsdaf|sqsubset|llcorner|vartheta|cupbrcap|lnapprox|Superset|SuchThat|succnsim|succneqq|angmsdag|biguplus|curlyvee|trpezium|Succeeds|NotTilde|bigwedge|angmsdah|angrtvbd|triminus|cwconint|fpartint|lrcorner|smeparsl|subseteq|urcorner|lurdshar|laemptyv|DDotrahd|approxeq|ldrushar|awconint|mapstoup|backcong|shortmid|triangle|geqslant|gesdotol|timesbar|circledR|circledS|setminus|multimap|naturals|scpolint|ncongdot|RightTee|boxminus|gnapprox|boxtimes|andslope|thicksim|angmsdaa|varsigma|cirfnint|rtriltri|angmsdab|rppolint|angmsdac|barwedge|drbkarow|clubsuit|thetasym|bsolhsub|capbrcup|dzigrarr|doteqdot|DotEqual|dotminus|UnderBar|NotEqual|realpart|otimesas|ulcorner|hksearow|hkswarow|parallel|PartialD|elinters|emptyset|plusacir|bbrktbrk|angmsdad|pointint|bigoplus|angmsdae|Precedes|bigsqcup|varkappa|notindot|supseteq|precneqq|precnsim|profalar|profline|profsurf|leqslant|lesdotor|raemptyv|subplus|notnivb|notnivc|subrarr|zigrarr|vzigzag|submult|subedot|Element|between|cirscir|larrbfs|larrsim|lotimes|lbrksld|lbrkslu|lozenge|ldrdhar|dbkarow|bigcirc|epsilon|simrarr|simplus|ltquest|Epsilon|luruhar|gtquest|maltese|npolint|eqcolon|npreceq|bigodot|ddagger|gtrless|bnequiv|harrcir|ddotseq|equivDD|backsim|demptyv|nsqsube|nsqsupe|Upsilon|nsubset|upsilon|minusdu|nsucceq|swarrow|nsupset|coloneq|searrow|boxplus|napprox|natural|asympeq|alefsym|congdot|nearrow|bigstar|diamond|supplus|tritime|LeftTee|nvinfin|triplus|NewLine|nvltrie|nvrtrie|nwarrow|nexists|Diamond|ruluhar|Implies|supmult|angzarr|suplarr|suphsub|questeq|because|digamma|Because|olcross|bemptyv|omicron|Omicron|rotimes|NoBreak|intprod|angrtvb|orderof|uwangle|suphsol|lesdoto|orslope|DownTee|realine|cudarrl|rdldhar|OverBar|supedot|lessdot|supdsub|topfork|succsim|rbrkslu|rbrksld|pertenk|cudarrr|isindot|planckh|lessgtr|pluscir|gesdoto|plussim|plustwo|lesssim|cularrp|rarrsim|Cayleys|notinva|notinvb|notinvc|UpArrow|Uparrow|uparrow|NotLess|dwangle|precsim|Product|curarrm|Cconint|dotplus|rarrbfs|ccupssm|Cedilla|cemptyv|notniva|quatint|frac35|frac38|frac45|frac56|frac58|frac78|tridot|xoplus|gacute|gammad|Gammad|lfisht|lfloor|bigcup|sqsupe|gbreve|Gbreve|lharul|sqsube|sqcups|Gcedil|apacir|llhard|lmidot|Lmidot|lmoust|andand|sqcaps|approx|Abreve|spades|circeq|tprime|divide|topcir|Assign|topbot|gesdot|divonx|xuplus|timesd|gesles|atilde|solbar|SOFTcy|loplus|timesb|lowast|lowbar|dlcorn|dlcrop|softcy|dollar|lparlt|thksim|lrhard|Atilde|lsaquo|smashp|bigvee|thinsp|wreath|bkarow|lsquor|lstrok|Lstrok|lthree|ltimes|ltlarr|DotDot|simdot|ltrPar|weierp|xsqcup|angmsd|sigmav|sigmaf|zeetrf|Zcaron|zcaron|mapsto|vsupne|thetav|cirmid|marker|mcomma|Zacute|vsubnE|there4|gtlPar|vsubne|bottom|gtrarr|SHCHcy|shchcy|midast|midcir|middot|minusb|minusd|gtrdot|bowtie|sfrown|mnplus|models|colone|seswar|Colone|mstpos|searhk|gtrsim|nacute|Nacute|boxbox|telrec|hairsp|Tcedil|nbumpe|scnsim|ncaron|Ncaron|ncedil|Ncedil|hamilt|Scedil|nearhk|hardcy|HARDcy|tcedil|Tcaron|commat|nequiv|nesear|tcaron|target|hearts|nexist|varrho|scedil|Scaron|scaron|hellip|Sacute|sacute|hercon|swnwar|compfn|rtimes|rthree|rsquor|rsaquo|zacute|wedgeq|homtht|barvee|barwed|Barwed|rpargt|horbar|conint|swarhk|roplus|nltrie|hslash|hstrok|Hstrok|rmoust|Conint|bprime|hybull|hyphen|iacute|Iacute|supsup|supsub|supsim|varphi|coprod|brvbar|agrave|Supset|supset|igrave|Igrave|notinE|Agrave|iiiint|iinfin|copysr|wedbar|Verbar|vangrt|becaus|incare|verbar|inodot|bullet|drcorn|intcal|drcrop|cularr|vellip|Utilde|bumpeq|cupcap|dstrok|Dstrok|CupCap|cupcup|cupdot|eacute|Eacute|supdot|iquest|easter|ecaron|Ecaron|ecolon|isinsv|utilde|itilde|Itilde|curarr|succeq|Bumpeq|cacute|ulcrop|nparsl|Cacute|nprcue|egrave|Egrave|nrarrc|nrarrw|subsup|subsub|nrtrie|jsercy|nsccue|Jsercy|kappav|kcedil|Kcedil|subsim|ulcorn|nsimeq|egsdot|veebar|kgreen|capand|elsdot|Subset|subset|curren|aacute|lacute|Lacute|emptyv|ntilde|Ntilde|lagran|lambda|Lambda|capcap|Ugrave|langle|subdot|emsp13|numero|emsp14|nvdash|nvDash|nVdash|nVDash|ugrave|ufisht|nvHarr|larrfs|nvlArr|larrhk|larrlp|larrpl|nvrArr|Udblac|nwarhk|larrtl|nwnear|oacute|Oacute|latail|lAtail|sstarf|lbrace|odblac|Odblac|lbrack|udblac|odsold|eparsl|lcaron|Lcaron|ograve|Ograve|lcedil|Lcedil|Aacute|ssmile|ssetmn|squarf|ldquor|capcup|ominus|cylcty|rharul|eqcirc|dagger|rfloor|rfisht|Dagger|daleth|equals|origof|capdot|equest|dcaron|Dcaron|rdquor|oslash|Oslash|otilde|Otilde|otimes|Otimes|urcrop|Ubreve|ubreve|Yacute|Uacute|uacute|Rcedil|rcedil|urcorn|parsim|Rcaron|Vdashl|rcaron|Tstrok|percnt|period|permil|Exists|yacute|rbrack|rbrace|phmmat|ccaron|Ccaron|planck|ccedil|plankv|tstrok|female|plusdo|plusdu|ffilig|plusmn|ffllig|Ccedil|rAtail|dfisht|bernou|ratail|Rarrtl|rarrtl|angsph|rarrpl|rarrlp|rarrhk|xwedge|xotime|forall|ForAll|Vvdash|vsupnE|preceq|bigcap|frac12|frac13|frac14|primes|rarrfs|prnsim|frac15|Square|frac16|square|lesdot|frac18|frac23|propto|prurel|rarrap|rangle|puncsp|frac25|Racute|qprime|racute|lesges|frac34|abreve|AElig|eqsim|utdot|setmn|urtri|Equal|Uring|seArr|uring|searr|dashv|Dashv|mumap|nabla|iogon|Iogon|sdote|sdotb|scsim|napid|napos|equiv|natur|Acirc|dblac|erarr|nbump|iprod|erDot|ucirc|awint|esdot|angrt|ncong|isinE|scnap|Scirc|scirc|ndash|isins|Ubrcy|nearr|neArr|isinv|nedot|ubrcy|acute|Ycirc|iukcy|Iukcy|xutri|nesim|caret|jcirc|Jcirc|caron|twixt|ddarr|sccue|exist|jmath|sbquo|ngeqq|angst|ccaps|lceil|ngsim|UpTee|delta|Delta|rtrif|nharr|nhArr|nhpar|rtrie|jukcy|Jukcy|kappa|rsquo|Kappa|nlarr|nlArr|TSHcy|rrarr|aogon|Aogon|fflig|xrarr|tshcy|ccirc|nleqq|filig|upsih|nless|dharl|nlsim|fjlig|ropar|nltri|dharr|robrk|roarr|fllig|fltns|roang|rnmid|subnE|subne|lAarr|trisb|Ccirc|acirc|ccups|blank|VDash|forkv|Vdash|langd|cedil|blk12|blk14|laquo|strns|diams|notin|vDash|larrb|blk34|block|disin|uplus|vdash|vBarv|aelig|starf|Wedge|check|xrArr|lates|lbarr|lBarr|notni|lbbrk|bcong|frasl|lbrke|frown|vrtri|vprop|vnsup|gamma|Gamma|wedge|xodot|bdquo|srarr|doteq|ldquo|boxdl|boxdL|gcirc|Gcirc|boxDl|boxDL|boxdr|boxdR|boxDr|TRADE|trade|rlhar|boxDR|vnsub|npart|vltri|rlarr|boxhd|boxhD|nprec|gescc|nrarr|nrArr|boxHd|boxHD|boxhu|boxhU|nrtri|boxHu|clubs|boxHU|times|colon|Colon|gimel|xlArr|Tilde|nsime|tilde|nsmid|nspar|THORN|thorn|xlarr|nsube|nsubE|thkap|xhArr|comma|nsucc|boxul|boxuL|nsupe|nsupE|gneqq|gnsim|boxUl|boxUL|grave|boxur|boxuR|boxUr|boxUR|lescc|angle|bepsi|boxvh|varpi|boxvH|numsp|Theta|gsime|gsiml|theta|boxVh|boxVH|boxvl|gtcir|gtdot|boxvL|boxVl|boxVL|crarr|cross|Cross|nvsim|boxvr|nwarr|nwArr|sqsup|dtdot|Uogon|lhard|lharu|dtrif|ocirc|Ocirc|lhblk|duarr|odash|sqsub|Hacek|sqcup|llarr|duhar|oelig|OElig|ofcir|boxvR|uogon|lltri|boxVr|csube|uuarr|ohbar|csupe|ctdot|olarr|olcir|harrw|oline|sqcap|omacr|Omacr|omega|Omega|boxVR|aleph|lneqq|lnsim|loang|loarr|rharu|lobrk|hcirc|operp|oplus|rhard|Hcirc|orarr|Union|order|ecirc|Ecirc|cuepr|szlig|cuesc|breve|reals|eDDot|Breve|hoarr|lopar|utrif|rdquo|Umacr|umacr|efDot|swArr|ultri|alpha|rceil|ovbar|swarr|Wcirc|wcirc|smtes|smile|bsemi|lrarr|aring|parsl|lrhar|bsime|uhblk|lrtri|cupor|Aring|uharr|uharl|slarr|rbrke|bsolb|lsime|rbbrk|RBarr|lsimg|phone|rBarr|rbarr|icirc|lsquo|Icirc|emacr|Emacr|ratio|simne|plusb|simlE|simgE|simeq|pluse|ltcir|ltdot|empty|xharr|xdtri|iexcl|Alpha|ltrie|rarrw|pound|ltrif|xcirc|bumpe|prcue|bumpE|asymp|amacr|cuvee|Sigma|sigma|iiint|udhar|iiota|ijlig|IJlig|supnE|imacr|Imacr|prime|Prime|image|prnap|eogon|Eogon|rarrc|mdash|mDDot|cuwed|imath|supne|imped|Amacr|udarr|prsim|micro|rarrb|cwint|raquo|infin|eplus|range|rangd|Ucirc|radic|minus|amalg|veeeq|rAarr|epsiv|ycirc|quest|sharp|quot|zwnj|Qscr|race|qscr|Qopf|qopf|qint|rang|Rang|Zscr|zscr|Zopf|zopf|rarr|rArr|Rarr|Pscr|pscr|prop|prod|prnE|prec|ZHcy|zhcy|prap|Zeta|zeta|Popf|popf|Zdot|plus|zdot|Yuml|yuml|phiv|YUcy|yucy|Yscr|yscr|perp|Yopf|yopf|part|para|YIcy|Ouml|rcub|yicy|YAcy|rdca|ouml|osol|Oscr|rdsh|yacy|real|oscr|xvee|andd|rect|andv|Xscr|oror|ordm|ordf|xscr|ange|aopf|Aopf|rHar|Xopf|opar|Oopf|xopf|xnis|rhov|oopf|omid|xmap|oint|apid|apos|ogon|ascr|Ascr|odot|odiv|xcup|xcap|ocir|oast|nvlt|nvle|nvgt|nvge|nvap|Wscr|wscr|auml|ntlg|ntgl|nsup|nsub|nsim|Nscr|nscr|nsce|Wopf|ring|npre|wopf|npar|Auml|Barv|bbrk|Nopf|nopf|nmid|nLtv|beta|ropf|Ropf|Beta|beth|nles|rpar|nleq|bnot|bNot|nldr|NJcy|rscr|Rscr|Vscr|vscr|rsqb|njcy|bopf|nisd|Bopf|rtri|Vopf|nGtv|ngtr|vopf|boxh|boxH|boxv|nges|ngeq|boxV|bscr|scap|Bscr|bsim|Vert|vert|bsol|bull|bump|caps|cdot|ncup|scnE|ncap|nbsp|napE|Cdot|cent|sdot|Vbar|nang|vBar|chcy|Mscr|mscr|sect|semi|CHcy|Mopf|mopf|sext|circ|cire|mldr|mlcp|cirE|comp|shcy|SHcy|vArr|varr|cong|copf|Copf|copy|COPY|malt|male|macr|lvnE|cscr|ltri|sime|ltcc|simg|Cscr|siml|csub|Uuml|lsqb|lsim|uuml|csup|Lscr|lscr|utri|smid|lpar|cups|smte|lozf|darr|Lopf|Uscr|solb|lopf|sopf|Sopf|lneq|uscr|spar|dArr|lnap|Darr|dash|Sqrt|LJcy|ljcy|lHar|dHar|Upsi|upsi|diam|lesg|djcy|DJcy|leqq|dopf|Dopf|dscr|Dscr|dscy|ldsh|ldca|squf|DScy|sscr|Sscr|dsol|lcub|late|star|Star|Uopf|Larr|lArr|larr|uopf|dtri|dzcy|sube|subE|Lang|lang|Kscr|kscr|Kopf|kopf|KJcy|kjcy|KHcy|khcy|DZcy|ecir|edot|eDot|Jscr|jscr|succ|Jopf|jopf|Edot|uHar|emsp|ensp|Iuml|iuml|eopf|isin|Iscr|iscr|Eopf|epar|sung|epsi|escr|sup1|sup2|sup3|Iota|iota|supe|supE|Iopf|iopf|IOcy|iocy|Escr|esim|Esim|imof|Uarr|QUOT|uArr|uarr|euml|IEcy|iecy|Idot|Euml|euro|excl|Hscr|hscr|Hopf|hopf|TScy|tscy|Tscr|hbar|tscr|flat|tbrk|fnof|hArr|harr|half|fopf|Fopf|tdot|gvnE|fork|trie|gtcc|fscr|Fscr|gdot|gsim|Gscr|gscr|Gopf|gopf|gneq|Gdot|tosa|gnap|Topf|topf|geqq|toea|GJcy|gjcy|tint|gesl|mid|Sfr|ggg|top|ges|gla|glE|glj|geq|gne|gEl|gel|gnE|Gcy|gcy|gap|Tfr|tfr|Tcy|tcy|Hat|Tau|Ffr|tau|Tab|hfr|Hfr|ffr|Fcy|fcy|icy|Icy|iff|ETH|eth|ifr|Ifr|Eta|eta|int|Int|Sup|sup|ucy|Ucy|Sum|sum|jcy|ENG|ufr|Ufr|eng|Jcy|jfr|els|ell|egs|Efr|efr|Jfr|uml|kcy|Kcy|Ecy|ecy|kfr|Kfr|lap|Sub|sub|lat|lcy|Lcy|leg|Dot|dot|lEg|leq|les|squ|div|die|lfr|Lfr|lgE|Dfr|dfr|Del|deg|Dcy|dcy|lne|lnE|sol|loz|smt|Cup|lrm|cup|lsh|Lsh|sim|shy|map|Map|mcy|Mcy|mfr|Mfr|mho|gfr|Gfr|sfr|cir|Chi|chi|nap|Cfr|vcy|Vcy|cfr|Scy|scy|ncy|Ncy|vee|Vee|Cap|cap|nfr|scE|sce|Nfr|nge|ngE|nGg|vfr|Vfr|ngt|bot|nGt|nis|niv|Rsh|rsh|nle|nlE|bne|Bfr|bfr|nLl|nlt|nLt|Bcy|bcy|not|Not|rlm|wfr|Wfr|npr|nsc|num|ocy|ast|Ocy|ofr|xfr|Xfr|Ofr|ogt|ohm|apE|olt|Rho|ape|rho|Rfr|rfr|ord|REG|ang|reg|orv|And|and|AMP|Rcy|amp|Afr|ycy|Ycy|yen|yfr|Yfr|rcy|par|pcy|Pcy|pfr|Pfr|phi|Phi|afr|Acy|acy|zcy|Zcy|piv|acE|acd|zfr|Zfr|pre|prE|psi|Psi|qfr|Qfr|zwj|Or|ge|Gg|gt|gg|el|oS|lt|Lt|LT|Re|lg|gl|eg|ne|Im|it|le|DD|wp|wr|nu|Nu|dd|lE|Sc|sc|pi|Pi|ee|af|ll|Ll|rx|gE|xi|pm|Xi|ic|pr|Pr|in|ni|mp|mu|ac|Mu|or|ap|Gt|GT|ii);|&(Aacute|Agrave|Atilde|Ccedil|Eacute|Egrave|Iacute|Igrave|Ntilde|Oacute|Ograve|Oslash|Otilde|Uacute|Ugrave|Yacute|aacute|agrave|atilde|brvbar|ccedil|curren|divide|eacute|egrave|frac12|frac14|frac34|iacute|igrave|iquest|middot|ntilde|oacute|ograve|oslash|otilde|plusmn|uacute|ugrave|yacute|AElig|Acirc|Aring|Ecirc|Icirc|Ocirc|THORN|Ucirc|acirc|acute|aelig|aring|cedil|ecirc|icirc|iexcl|laquo|micro|ocirc|pound|raquo|szlig|thorn|times|ucirc|Auml|COPY|Euml|Iuml|Ouml|QUOT|Uuml|auml|cent|copy|euml|iuml|macr|nbsp|ordf|ordm|ouml|para|quot|sect|sup1|sup2|sup3|uuml|yuml|AMP|ETH|REG|amp|deg|eth|not|reg|shy|uml|yen|GT|LT|gt|lt)(?!;)([=a-zA-Z0-9]?)|&#([0-9]+)(;?)|&#[xX]([a-fA-F0-9]+)(;?)|&([0-9a-zA-Z]+)/g;
		var decodeMap = {'aacute':'\xE1','Aacute':'\xC1','abreve':'\u0103','Abreve':'\u0102','ac':'\u223E','acd':'\u223F','acE':'\u223E\u0333','acirc':'\xE2','Acirc':'\xC2','acute':'\xB4','acy':'\u0430','Acy':'\u0410','aelig':'\xE6','AElig':'\xC6','af':'\u2061','afr':'\uD835\uDD1E','Afr':'\uD835\uDD04','agrave':'\xE0','Agrave':'\xC0','alefsym':'\u2135','aleph':'\u2135','alpha':'\u03B1','Alpha':'\u0391','amacr':'\u0101','Amacr':'\u0100','amalg':'\u2A3F','amp':'&','AMP':'&','and':'\u2227','And':'\u2A53','andand':'\u2A55','andd':'\u2A5C','andslope':'\u2A58','andv':'\u2A5A','ang':'\u2220','ange':'\u29A4','angle':'\u2220','angmsd':'\u2221','angmsdaa':'\u29A8','angmsdab':'\u29A9','angmsdac':'\u29AA','angmsdad':'\u29AB','angmsdae':'\u29AC','angmsdaf':'\u29AD','angmsdag':'\u29AE','angmsdah':'\u29AF','angrt':'\u221F','angrtvb':'\u22BE','angrtvbd':'\u299D','angsph':'\u2222','angst':'\xC5','angzarr':'\u237C','aogon':'\u0105','Aogon':'\u0104','aopf':'\uD835\uDD52','Aopf':'\uD835\uDD38','ap':'\u2248','apacir':'\u2A6F','ape':'\u224A','apE':'\u2A70','apid':'\u224B','apos':'\'','ApplyFunction':'\u2061','approx':'\u2248','approxeq':'\u224A','aring':'\xE5','Aring':'\xC5','ascr':'\uD835\uDCB6','Ascr':'\uD835\uDC9C','Assign':'\u2254','ast':'*','asymp':'\u2248','asympeq':'\u224D','atilde':'\xE3','Atilde':'\xC3','auml':'\xE4','Auml':'\xC4','awconint':'\u2233','awint':'\u2A11','backcong':'\u224C','backepsilon':'\u03F6','backprime':'\u2035','backsim':'\u223D','backsimeq':'\u22CD','Backslash':'\u2216','Barv':'\u2AE7','barvee':'\u22BD','barwed':'\u2305','Barwed':'\u2306','barwedge':'\u2305','bbrk':'\u23B5','bbrktbrk':'\u23B6','bcong':'\u224C','bcy':'\u0431','Bcy':'\u0411','bdquo':'\u201E','becaus':'\u2235','because':'\u2235','Because':'\u2235','bemptyv':'\u29B0','bepsi':'\u03F6','bernou':'\u212C','Bernoullis':'\u212C','beta':'\u03B2','Beta':'\u0392','beth':'\u2136','between':'\u226C','bfr':'\uD835\uDD1F','Bfr':'\uD835\uDD05','bigcap':'\u22C2','bigcirc':'\u25EF','bigcup':'\u22C3','bigodot':'\u2A00','bigoplus':'\u2A01','bigotimes':'\u2A02','bigsqcup':'\u2A06','bigstar':'\u2605','bigtriangledown':'\u25BD','bigtriangleup':'\u25B3','biguplus':'\u2A04','bigvee':'\u22C1','bigwedge':'\u22C0','bkarow':'\u290D','blacklozenge':'\u29EB','blacksquare':'\u25AA','blacktriangle':'\u25B4','blacktriangledown':'\u25BE','blacktriangleleft':'\u25C2','blacktriangleright':'\u25B8','blank':'\u2423','blk12':'\u2592','blk14':'\u2591','blk34':'\u2593','block':'\u2588','bne':'=\u20E5','bnequiv':'\u2261\u20E5','bnot':'\u2310','bNot':'\u2AED','bopf':'\uD835\uDD53','Bopf':'\uD835\uDD39','bot':'\u22A5','bottom':'\u22A5','bowtie':'\u22C8','boxbox':'\u29C9','boxdl':'\u2510','boxdL':'\u2555','boxDl':'\u2556','boxDL':'\u2557','boxdr':'\u250C','boxdR':'\u2552','boxDr':'\u2553','boxDR':'\u2554','boxh':'\u2500','boxH':'\u2550','boxhd':'\u252C','boxhD':'\u2565','boxHd':'\u2564','boxHD':'\u2566','boxhu':'\u2534','boxhU':'\u2568','boxHu':'\u2567','boxHU':'\u2569','boxminus':'\u229F','boxplus':'\u229E','boxtimes':'\u22A0','boxul':'\u2518','boxuL':'\u255B','boxUl':'\u255C','boxUL':'\u255D','boxur':'\u2514','boxuR':'\u2558','boxUr':'\u2559','boxUR':'\u255A','boxv':'\u2502','boxV':'\u2551','boxvh':'\u253C','boxvH':'\u256A','boxVh':'\u256B','boxVH':'\u256C','boxvl':'\u2524','boxvL':'\u2561','boxVl':'\u2562','boxVL':'\u2563','boxvr':'\u251C','boxvR':'\u255E','boxVr':'\u255F','boxVR':'\u2560','bprime':'\u2035','breve':'\u02D8','Breve':'\u02D8','brvbar':'\xA6','bscr':'\uD835\uDCB7','Bscr':'\u212C','bsemi':'\u204F','bsim':'\u223D','bsime':'\u22CD','bsol':'\\','bsolb':'\u29C5','bsolhsub':'\u27C8','bull':'\u2022','bullet':'\u2022','bump':'\u224E','bumpe':'\u224F','bumpE':'\u2AAE','bumpeq':'\u224F','Bumpeq':'\u224E','cacute':'\u0107','Cacute':'\u0106','cap':'\u2229','Cap':'\u22D2','capand':'\u2A44','capbrcup':'\u2A49','capcap':'\u2A4B','capcup':'\u2A47','capdot':'\u2A40','CapitalDifferentialD':'\u2145','caps':'\u2229\uFE00','caret':'\u2041','caron':'\u02C7','Cayleys':'\u212D','ccaps':'\u2A4D','ccaron':'\u010D','Ccaron':'\u010C','ccedil':'\xE7','Ccedil':'\xC7','ccirc':'\u0109','Ccirc':'\u0108','Cconint':'\u2230','ccups':'\u2A4C','ccupssm':'\u2A50','cdot':'\u010B','Cdot':'\u010A','cedil':'\xB8','Cedilla':'\xB8','cemptyv':'\u29B2','cent':'\xA2','centerdot':'\xB7','CenterDot':'\xB7','cfr':'\uD835\uDD20','Cfr':'\u212D','chcy':'\u0447','CHcy':'\u0427','check':'\u2713','checkmark':'\u2713','chi':'\u03C7','Chi':'\u03A7','cir':'\u25CB','circ':'\u02C6','circeq':'\u2257','circlearrowleft':'\u21BA','circlearrowright':'\u21BB','circledast':'\u229B','circledcirc':'\u229A','circleddash':'\u229D','CircleDot':'\u2299','circledR':'\xAE','circledS':'\u24C8','CircleMinus':'\u2296','CirclePlus':'\u2295','CircleTimes':'\u2297','cire':'\u2257','cirE':'\u29C3','cirfnint':'\u2A10','cirmid':'\u2AEF','cirscir':'\u29C2','ClockwiseContourIntegral':'\u2232','CloseCurlyDoubleQuote':'\u201D','CloseCurlyQuote':'\u2019','clubs':'\u2663','clubsuit':'\u2663','colon':':','Colon':'\u2237','colone':'\u2254','Colone':'\u2A74','coloneq':'\u2254','comma':',','commat':'@','comp':'\u2201','compfn':'\u2218','complement':'\u2201','complexes':'\u2102','cong':'\u2245','congdot':'\u2A6D','Congruent':'\u2261','conint':'\u222E','Conint':'\u222F','ContourIntegral':'\u222E','copf':'\uD835\uDD54','Copf':'\u2102','coprod':'\u2210','Coproduct':'\u2210','copy':'\xA9','COPY':'\xA9','copysr':'\u2117','CounterClockwiseContourIntegral':'\u2233','crarr':'\u21B5','cross':'\u2717','Cross':'\u2A2F','cscr':'\uD835\uDCB8','Cscr':'\uD835\uDC9E','csub':'\u2ACF','csube':'\u2AD1','csup':'\u2AD0','csupe':'\u2AD2','ctdot':'\u22EF','cudarrl':'\u2938','cudarrr':'\u2935','cuepr':'\u22DE','cuesc':'\u22DF','cularr':'\u21B6','cularrp':'\u293D','cup':'\u222A','Cup':'\u22D3','cupbrcap':'\u2A48','cupcap':'\u2A46','CupCap':'\u224D','cupcup':'\u2A4A','cupdot':'\u228D','cupor':'\u2A45','cups':'\u222A\uFE00','curarr':'\u21B7','curarrm':'\u293C','curlyeqprec':'\u22DE','curlyeqsucc':'\u22DF','curlyvee':'\u22CE','curlywedge':'\u22CF','curren':'\xA4','curvearrowleft':'\u21B6','curvearrowright':'\u21B7','cuvee':'\u22CE','cuwed':'\u22CF','cwconint':'\u2232','cwint':'\u2231','cylcty':'\u232D','dagger':'\u2020','Dagger':'\u2021','daleth':'\u2138','darr':'\u2193','dArr':'\u21D3','Darr':'\u21A1','dash':'\u2010','dashv':'\u22A3','Dashv':'\u2AE4','dbkarow':'\u290F','dblac':'\u02DD','dcaron':'\u010F','Dcaron':'\u010E','dcy':'\u0434','Dcy':'\u0414','dd':'\u2146','DD':'\u2145','ddagger':'\u2021','ddarr':'\u21CA','DDotrahd':'\u2911','ddotseq':'\u2A77','deg':'\xB0','Del':'\u2207','delta':'\u03B4','Delta':'\u0394','demptyv':'\u29B1','dfisht':'\u297F','dfr':'\uD835\uDD21','Dfr':'\uD835\uDD07','dHar':'\u2965','dharl':'\u21C3','dharr':'\u21C2','DiacriticalAcute':'\xB4','DiacriticalDot':'\u02D9','DiacriticalDoubleAcute':'\u02DD','DiacriticalGrave':'`','DiacriticalTilde':'\u02DC','diam':'\u22C4','diamond':'\u22C4','Diamond':'\u22C4','diamondsuit':'\u2666','diams':'\u2666','die':'\xA8','DifferentialD':'\u2146','digamma':'\u03DD','disin':'\u22F2','div':'\xF7','divide':'\xF7','divideontimes':'\u22C7','divonx':'\u22C7','djcy':'\u0452','DJcy':'\u0402','dlcorn':'\u231E','dlcrop':'\u230D','dollar':'$','dopf':'\uD835\uDD55','Dopf':'\uD835\uDD3B','dot':'\u02D9','Dot':'\xA8','DotDot':'\u20DC','doteq':'\u2250','doteqdot':'\u2251','DotEqual':'\u2250','dotminus':'\u2238','dotplus':'\u2214','dotsquare':'\u22A1','doublebarwedge':'\u2306','DoubleContourIntegral':'\u222F','DoubleDot':'\xA8','DoubleDownArrow':'\u21D3','DoubleLeftArrow':'\u21D0','DoubleLeftRightArrow':'\u21D4','DoubleLeftTee':'\u2AE4','DoubleLongLeftArrow':'\u27F8','DoubleLongLeftRightArrow':'\u27FA','DoubleLongRightArrow':'\u27F9','DoubleRightArrow':'\u21D2','DoubleRightTee':'\u22A8','DoubleUpArrow':'\u21D1','DoubleUpDownArrow':'\u21D5','DoubleVerticalBar':'\u2225','downarrow':'\u2193','Downarrow':'\u21D3','DownArrow':'\u2193','DownArrowBar':'\u2913','DownArrowUpArrow':'\u21F5','DownBreve':'\u0311','downdownarrows':'\u21CA','downharpoonleft':'\u21C3','downharpoonright':'\u21C2','DownLeftRightVector':'\u2950','DownLeftTeeVector':'\u295E','DownLeftVector':'\u21BD','DownLeftVectorBar':'\u2956','DownRightTeeVector':'\u295F','DownRightVector':'\u21C1','DownRightVectorBar':'\u2957','DownTee':'\u22A4','DownTeeArrow':'\u21A7','drbkarow':'\u2910','drcorn':'\u231F','drcrop':'\u230C','dscr':'\uD835\uDCB9','Dscr':'\uD835\uDC9F','dscy':'\u0455','DScy':'\u0405','dsol':'\u29F6','dstrok':'\u0111','Dstrok':'\u0110','dtdot':'\u22F1','dtri':'\u25BF','dtrif':'\u25BE','duarr':'\u21F5','duhar':'\u296F','dwangle':'\u29A6','dzcy':'\u045F','DZcy':'\u040F','dzigrarr':'\u27FF','eacute':'\xE9','Eacute':'\xC9','easter':'\u2A6E','ecaron':'\u011B','Ecaron':'\u011A','ecir':'\u2256','ecirc':'\xEA','Ecirc':'\xCA','ecolon':'\u2255','ecy':'\u044D','Ecy':'\u042D','eDDot':'\u2A77','edot':'\u0117','eDot':'\u2251','Edot':'\u0116','ee':'\u2147','efDot':'\u2252','efr':'\uD835\uDD22','Efr':'\uD835\uDD08','eg':'\u2A9A','egrave':'\xE8','Egrave':'\xC8','egs':'\u2A96','egsdot':'\u2A98','el':'\u2A99','Element':'\u2208','elinters':'\u23E7','ell':'\u2113','els':'\u2A95','elsdot':'\u2A97','emacr':'\u0113','Emacr':'\u0112','empty':'\u2205','emptyset':'\u2205','EmptySmallSquare':'\u25FB','emptyv':'\u2205','EmptyVerySmallSquare':'\u25AB','emsp':'\u2003','emsp13':'\u2004','emsp14':'\u2005','eng':'\u014B','ENG':'\u014A','ensp':'\u2002','eogon':'\u0119','Eogon':'\u0118','eopf':'\uD835\uDD56','Eopf':'\uD835\uDD3C','epar':'\u22D5','eparsl':'\u29E3','eplus':'\u2A71','epsi':'\u03B5','epsilon':'\u03B5','Epsilon':'\u0395','epsiv':'\u03F5','eqcirc':'\u2256','eqcolon':'\u2255','eqsim':'\u2242','eqslantgtr':'\u2A96','eqslantless':'\u2A95','Equal':'\u2A75','equals':'=','EqualTilde':'\u2242','equest':'\u225F','Equilibrium':'\u21CC','equiv':'\u2261','equivDD':'\u2A78','eqvparsl':'\u29E5','erarr':'\u2971','erDot':'\u2253','escr':'\u212F','Escr':'\u2130','esdot':'\u2250','esim':'\u2242','Esim':'\u2A73','eta':'\u03B7','Eta':'\u0397','eth':'\xF0','ETH':'\xD0','euml':'\xEB','Euml':'\xCB','euro':'\u20AC','excl':'!','exist':'\u2203','Exists':'\u2203','expectation':'\u2130','exponentiale':'\u2147','ExponentialE':'\u2147','fallingdotseq':'\u2252','fcy':'\u0444','Fcy':'\u0424','female':'\u2640','ffilig':'\uFB03','fflig':'\uFB00','ffllig':'\uFB04','ffr':'\uD835\uDD23','Ffr':'\uD835\uDD09','filig':'\uFB01','FilledSmallSquare':'\u25FC','FilledVerySmallSquare':'\u25AA','fjlig':'fj','flat':'\u266D','fllig':'\uFB02','fltns':'\u25B1','fnof':'\u0192','fopf':'\uD835\uDD57','Fopf':'\uD835\uDD3D','forall':'\u2200','ForAll':'\u2200','fork':'\u22D4','forkv':'\u2AD9','Fouriertrf':'\u2131','fpartint':'\u2A0D','frac12':'\xBD','frac13':'\u2153','frac14':'\xBC','frac15':'\u2155','frac16':'\u2159','frac18':'\u215B','frac23':'\u2154','frac25':'\u2156','frac34':'\xBE','frac35':'\u2157','frac38':'\u215C','frac45':'\u2158','frac56':'\u215A','frac58':'\u215D','frac78':'\u215E','frasl':'\u2044','frown':'\u2322','fscr':'\uD835\uDCBB','Fscr':'\u2131','gacute':'\u01F5','gamma':'\u03B3','Gamma':'\u0393','gammad':'\u03DD','Gammad':'\u03DC','gap':'\u2A86','gbreve':'\u011F','Gbreve':'\u011E','Gcedil':'\u0122','gcirc':'\u011D','Gcirc':'\u011C','gcy':'\u0433','Gcy':'\u0413','gdot':'\u0121','Gdot':'\u0120','ge':'\u2265','gE':'\u2267','gel':'\u22DB','gEl':'\u2A8C','geq':'\u2265','geqq':'\u2267','geqslant':'\u2A7E','ges':'\u2A7E','gescc':'\u2AA9','gesdot':'\u2A80','gesdoto':'\u2A82','gesdotol':'\u2A84','gesl':'\u22DB\uFE00','gesles':'\u2A94','gfr':'\uD835\uDD24','Gfr':'\uD835\uDD0A','gg':'\u226B','Gg':'\u22D9','ggg':'\u22D9','gimel':'\u2137','gjcy':'\u0453','GJcy':'\u0403','gl':'\u2277','gla':'\u2AA5','glE':'\u2A92','glj':'\u2AA4','gnap':'\u2A8A','gnapprox':'\u2A8A','gne':'\u2A88','gnE':'\u2269','gneq':'\u2A88','gneqq':'\u2269','gnsim':'\u22E7','gopf':'\uD835\uDD58','Gopf':'\uD835\uDD3E','grave':'`','GreaterEqual':'\u2265','GreaterEqualLess':'\u22DB','GreaterFullEqual':'\u2267','GreaterGreater':'\u2AA2','GreaterLess':'\u2277','GreaterSlantEqual':'\u2A7E','GreaterTilde':'\u2273','gscr':'\u210A','Gscr':'\uD835\uDCA2','gsim':'\u2273','gsime':'\u2A8E','gsiml':'\u2A90','gt':'>','Gt':'\u226B','GT':'>','gtcc':'\u2AA7','gtcir':'\u2A7A','gtdot':'\u22D7','gtlPar':'\u2995','gtquest':'\u2A7C','gtrapprox':'\u2A86','gtrarr':'\u2978','gtrdot':'\u22D7','gtreqless':'\u22DB','gtreqqless':'\u2A8C','gtrless':'\u2277','gtrsim':'\u2273','gvertneqq':'\u2269\uFE00','gvnE':'\u2269\uFE00','Hacek':'\u02C7','hairsp':'\u200A','half':'\xBD','hamilt':'\u210B','hardcy':'\u044A','HARDcy':'\u042A','harr':'\u2194','hArr':'\u21D4','harrcir':'\u2948','harrw':'\u21AD','Hat':'^','hbar':'\u210F','hcirc':'\u0125','Hcirc':'\u0124','hearts':'\u2665','heartsuit':'\u2665','hellip':'\u2026','hercon':'\u22B9','hfr':'\uD835\uDD25','Hfr':'\u210C','HilbertSpace':'\u210B','hksearow':'\u2925','hkswarow':'\u2926','hoarr':'\u21FF','homtht':'\u223B','hookleftarrow':'\u21A9','hookrightarrow':'\u21AA','hopf':'\uD835\uDD59','Hopf':'\u210D','horbar':'\u2015','HorizontalLine':'\u2500','hscr':'\uD835\uDCBD','Hscr':'\u210B','hslash':'\u210F','hstrok':'\u0127','Hstrok':'\u0126','HumpDownHump':'\u224E','HumpEqual':'\u224F','hybull':'\u2043','hyphen':'\u2010','iacute':'\xED','Iacute':'\xCD','ic':'\u2063','icirc':'\xEE','Icirc':'\xCE','icy':'\u0438','Icy':'\u0418','Idot':'\u0130','iecy':'\u0435','IEcy':'\u0415','iexcl':'\xA1','iff':'\u21D4','ifr':'\uD835\uDD26','Ifr':'\u2111','igrave':'\xEC','Igrave':'\xCC','ii':'\u2148','iiiint':'\u2A0C','iiint':'\u222D','iinfin':'\u29DC','iiota':'\u2129','ijlig':'\u0133','IJlig':'\u0132','Im':'\u2111','imacr':'\u012B','Imacr':'\u012A','image':'\u2111','ImaginaryI':'\u2148','imagline':'\u2110','imagpart':'\u2111','imath':'\u0131','imof':'\u22B7','imped':'\u01B5','Implies':'\u21D2','in':'\u2208','incare':'\u2105','infin':'\u221E','infintie':'\u29DD','inodot':'\u0131','int':'\u222B','Int':'\u222C','intcal':'\u22BA','integers':'\u2124','Integral':'\u222B','intercal':'\u22BA','Intersection':'\u22C2','intlarhk':'\u2A17','intprod':'\u2A3C','InvisibleComma':'\u2063','InvisibleTimes':'\u2062','iocy':'\u0451','IOcy':'\u0401','iogon':'\u012F','Iogon':'\u012E','iopf':'\uD835\uDD5A','Iopf':'\uD835\uDD40','iota':'\u03B9','Iota':'\u0399','iprod':'\u2A3C','iquest':'\xBF','iscr':'\uD835\uDCBE','Iscr':'\u2110','isin':'\u2208','isindot':'\u22F5','isinE':'\u22F9','isins':'\u22F4','isinsv':'\u22F3','isinv':'\u2208','it':'\u2062','itilde':'\u0129','Itilde':'\u0128','iukcy':'\u0456','Iukcy':'\u0406','iuml':'\xEF','Iuml':'\xCF','jcirc':'\u0135','Jcirc':'\u0134','jcy':'\u0439','Jcy':'\u0419','jfr':'\uD835\uDD27','Jfr':'\uD835\uDD0D','jmath':'\u0237','jopf':'\uD835\uDD5B','Jopf':'\uD835\uDD41','jscr':'\uD835\uDCBF','Jscr':'\uD835\uDCA5','jsercy':'\u0458','Jsercy':'\u0408','jukcy':'\u0454','Jukcy':'\u0404','kappa':'\u03BA','Kappa':'\u039A','kappav':'\u03F0','kcedil':'\u0137','Kcedil':'\u0136','kcy':'\u043A','Kcy':'\u041A','kfr':'\uD835\uDD28','Kfr':'\uD835\uDD0E','kgreen':'\u0138','khcy':'\u0445','KHcy':'\u0425','kjcy':'\u045C','KJcy':'\u040C','kopf':'\uD835\uDD5C','Kopf':'\uD835\uDD42','kscr':'\uD835\uDCC0','Kscr':'\uD835\uDCA6','lAarr':'\u21DA','lacute':'\u013A','Lacute':'\u0139','laemptyv':'\u29B4','lagran':'\u2112','lambda':'\u03BB','Lambda':'\u039B','lang':'\u27E8','Lang':'\u27EA','langd':'\u2991','langle':'\u27E8','lap':'\u2A85','Laplacetrf':'\u2112','laquo':'\xAB','larr':'\u2190','lArr':'\u21D0','Larr':'\u219E','larrb':'\u21E4','larrbfs':'\u291F','larrfs':'\u291D','larrhk':'\u21A9','larrlp':'\u21AB','larrpl':'\u2939','larrsim':'\u2973','larrtl':'\u21A2','lat':'\u2AAB','latail':'\u2919','lAtail':'\u291B','late':'\u2AAD','lates':'\u2AAD\uFE00','lbarr':'\u290C','lBarr':'\u290E','lbbrk':'\u2772','lbrace':'{','lbrack':'[','lbrke':'\u298B','lbrksld':'\u298F','lbrkslu':'\u298D','lcaron':'\u013E','Lcaron':'\u013D','lcedil':'\u013C','Lcedil':'\u013B','lceil':'\u2308','lcub':'{','lcy':'\u043B','Lcy':'\u041B','ldca':'\u2936','ldquo':'\u201C','ldquor':'\u201E','ldrdhar':'\u2967','ldrushar':'\u294B','ldsh':'\u21B2','le':'\u2264','lE':'\u2266','LeftAngleBracket':'\u27E8','leftarrow':'\u2190','Leftarrow':'\u21D0','LeftArrow':'\u2190','LeftArrowBar':'\u21E4','LeftArrowRightArrow':'\u21C6','leftarrowtail':'\u21A2','LeftCeiling':'\u2308','LeftDoubleBracket':'\u27E6','LeftDownTeeVector':'\u2961','LeftDownVector':'\u21C3','LeftDownVectorBar':'\u2959','LeftFloor':'\u230A','leftharpoondown':'\u21BD','leftharpoonup':'\u21BC','leftleftarrows':'\u21C7','leftrightarrow':'\u2194','Leftrightarrow':'\u21D4','LeftRightArrow':'\u2194','leftrightarrows':'\u21C6','leftrightharpoons':'\u21CB','leftrightsquigarrow':'\u21AD','LeftRightVector':'\u294E','LeftTee':'\u22A3','LeftTeeArrow':'\u21A4','LeftTeeVector':'\u295A','leftthreetimes':'\u22CB','LeftTriangle':'\u22B2','LeftTriangleBar':'\u29CF','LeftTriangleEqual':'\u22B4','LeftUpDownVector':'\u2951','LeftUpTeeVector':'\u2960','LeftUpVector':'\u21BF','LeftUpVectorBar':'\u2958','LeftVector':'\u21BC','LeftVectorBar':'\u2952','leg':'\u22DA','lEg':'\u2A8B','leq':'\u2264','leqq':'\u2266','leqslant':'\u2A7D','les':'\u2A7D','lescc':'\u2AA8','lesdot':'\u2A7F','lesdoto':'\u2A81','lesdotor':'\u2A83','lesg':'\u22DA\uFE00','lesges':'\u2A93','lessapprox':'\u2A85','lessdot':'\u22D6','lesseqgtr':'\u22DA','lesseqqgtr':'\u2A8B','LessEqualGreater':'\u22DA','LessFullEqual':'\u2266','LessGreater':'\u2276','lessgtr':'\u2276','LessLess':'\u2AA1','lesssim':'\u2272','LessSlantEqual':'\u2A7D','LessTilde':'\u2272','lfisht':'\u297C','lfloor':'\u230A','lfr':'\uD835\uDD29','Lfr':'\uD835\uDD0F','lg':'\u2276','lgE':'\u2A91','lHar':'\u2962','lhard':'\u21BD','lharu':'\u21BC','lharul':'\u296A','lhblk':'\u2584','ljcy':'\u0459','LJcy':'\u0409','ll':'\u226A','Ll':'\u22D8','llarr':'\u21C7','llcorner':'\u231E','Lleftarrow':'\u21DA','llhard':'\u296B','lltri':'\u25FA','lmidot':'\u0140','Lmidot':'\u013F','lmoust':'\u23B0','lmoustache':'\u23B0','lnap':'\u2A89','lnapprox':'\u2A89','lne':'\u2A87','lnE':'\u2268','lneq':'\u2A87','lneqq':'\u2268','lnsim':'\u22E6','loang':'\u27EC','loarr':'\u21FD','lobrk':'\u27E6','longleftarrow':'\u27F5','Longleftarrow':'\u27F8','LongLeftArrow':'\u27F5','longleftrightarrow':'\u27F7','Longleftrightarrow':'\u27FA','LongLeftRightArrow':'\u27F7','longmapsto':'\u27FC','longrightarrow':'\u27F6','Longrightarrow':'\u27F9','LongRightArrow':'\u27F6','looparrowleft':'\u21AB','looparrowright':'\u21AC','lopar':'\u2985','lopf':'\uD835\uDD5D','Lopf':'\uD835\uDD43','loplus':'\u2A2D','lotimes':'\u2A34','lowast':'\u2217','lowbar':'_','LowerLeftArrow':'\u2199','LowerRightArrow':'\u2198','loz':'\u25CA','lozenge':'\u25CA','lozf':'\u29EB','lpar':'(','lparlt':'\u2993','lrarr':'\u21C6','lrcorner':'\u231F','lrhar':'\u21CB','lrhard':'\u296D','lrm':'\u200E','lrtri':'\u22BF','lsaquo':'\u2039','lscr':'\uD835\uDCC1','Lscr':'\u2112','lsh':'\u21B0','Lsh':'\u21B0','lsim':'\u2272','lsime':'\u2A8D','lsimg':'\u2A8F','lsqb':'[','lsquo':'\u2018','lsquor':'\u201A','lstrok':'\u0142','Lstrok':'\u0141','lt':'<','Lt':'\u226A','LT':'<','ltcc':'\u2AA6','ltcir':'\u2A79','ltdot':'\u22D6','lthree':'\u22CB','ltimes':'\u22C9','ltlarr':'\u2976','ltquest':'\u2A7B','ltri':'\u25C3','ltrie':'\u22B4','ltrif':'\u25C2','ltrPar':'\u2996','lurdshar':'\u294A','luruhar':'\u2966','lvertneqq':'\u2268\uFE00','lvnE':'\u2268\uFE00','macr':'\xAF','male':'\u2642','malt':'\u2720','maltese':'\u2720','map':'\u21A6','Map':'\u2905','mapsto':'\u21A6','mapstodown':'\u21A7','mapstoleft':'\u21A4','mapstoup':'\u21A5','marker':'\u25AE','mcomma':'\u2A29','mcy':'\u043C','Mcy':'\u041C','mdash':'\u2014','mDDot':'\u223A','measuredangle':'\u2221','MediumSpace':'\u205F','Mellintrf':'\u2133','mfr':'\uD835\uDD2A','Mfr':'\uD835\uDD10','mho':'\u2127','micro':'\xB5','mid':'\u2223','midast':'*','midcir':'\u2AF0','middot':'\xB7','minus':'\u2212','minusb':'\u229F','minusd':'\u2238','minusdu':'\u2A2A','MinusPlus':'\u2213','mlcp':'\u2ADB','mldr':'\u2026','mnplus':'\u2213','models':'\u22A7','mopf':'\uD835\uDD5E','Mopf':'\uD835\uDD44','mp':'\u2213','mscr':'\uD835\uDCC2','Mscr':'\u2133','mstpos':'\u223E','mu':'\u03BC','Mu':'\u039C','multimap':'\u22B8','mumap':'\u22B8','nabla':'\u2207','nacute':'\u0144','Nacute':'\u0143','nang':'\u2220\u20D2','nap':'\u2249','napE':'\u2A70\u0338','napid':'\u224B\u0338','napos':'\u0149','napprox':'\u2249','natur':'\u266E','natural':'\u266E','naturals':'\u2115','nbsp':'\xA0','nbump':'\u224E\u0338','nbumpe':'\u224F\u0338','ncap':'\u2A43','ncaron':'\u0148','Ncaron':'\u0147','ncedil':'\u0146','Ncedil':'\u0145','ncong':'\u2247','ncongdot':'\u2A6D\u0338','ncup':'\u2A42','ncy':'\u043D','Ncy':'\u041D','ndash':'\u2013','ne':'\u2260','nearhk':'\u2924','nearr':'\u2197','neArr':'\u21D7','nearrow':'\u2197','nedot':'\u2250\u0338','NegativeMediumSpace':'\u200B','NegativeThickSpace':'\u200B','NegativeThinSpace':'\u200B','NegativeVeryThinSpace':'\u200B','nequiv':'\u2262','nesear':'\u2928','nesim':'\u2242\u0338','NestedGreaterGreater':'\u226B','NestedLessLess':'\u226A','NewLine':'\n','nexist':'\u2204','nexists':'\u2204','nfr':'\uD835\uDD2B','Nfr':'\uD835\uDD11','nge':'\u2271','ngE':'\u2267\u0338','ngeq':'\u2271','ngeqq':'\u2267\u0338','ngeqslant':'\u2A7E\u0338','nges':'\u2A7E\u0338','nGg':'\u22D9\u0338','ngsim':'\u2275','ngt':'\u226F','nGt':'\u226B\u20D2','ngtr':'\u226F','nGtv':'\u226B\u0338','nharr':'\u21AE','nhArr':'\u21CE','nhpar':'\u2AF2','ni':'\u220B','nis':'\u22FC','nisd':'\u22FA','niv':'\u220B','njcy':'\u045A','NJcy':'\u040A','nlarr':'\u219A','nlArr':'\u21CD','nldr':'\u2025','nle':'\u2270','nlE':'\u2266\u0338','nleftarrow':'\u219A','nLeftarrow':'\u21CD','nleftrightarrow':'\u21AE','nLeftrightarrow':'\u21CE','nleq':'\u2270','nleqq':'\u2266\u0338','nleqslant':'\u2A7D\u0338','nles':'\u2A7D\u0338','nless':'\u226E','nLl':'\u22D8\u0338','nlsim':'\u2274','nlt':'\u226E','nLt':'\u226A\u20D2','nltri':'\u22EA','nltrie':'\u22EC','nLtv':'\u226A\u0338','nmid':'\u2224','NoBreak':'\u2060','NonBreakingSpace':'\xA0','nopf':'\uD835\uDD5F','Nopf':'\u2115','not':'\xAC','Not':'\u2AEC','NotCongruent':'\u2262','NotCupCap':'\u226D','NotDoubleVerticalBar':'\u2226','NotElement':'\u2209','NotEqual':'\u2260','NotEqualTilde':'\u2242\u0338','NotExists':'\u2204','NotGreater':'\u226F','NotGreaterEqual':'\u2271','NotGreaterFullEqual':'\u2267\u0338','NotGreaterGreater':'\u226B\u0338','NotGreaterLess':'\u2279','NotGreaterSlantEqual':'\u2A7E\u0338','NotGreaterTilde':'\u2275','NotHumpDownHump':'\u224E\u0338','NotHumpEqual':'\u224F\u0338','notin':'\u2209','notindot':'\u22F5\u0338','notinE':'\u22F9\u0338','notinva':'\u2209','notinvb':'\u22F7','notinvc':'\u22F6','NotLeftTriangle':'\u22EA','NotLeftTriangleBar':'\u29CF\u0338','NotLeftTriangleEqual':'\u22EC','NotLess':'\u226E','NotLessEqual':'\u2270','NotLessGreater':'\u2278','NotLessLess':'\u226A\u0338','NotLessSlantEqual':'\u2A7D\u0338','NotLessTilde':'\u2274','NotNestedGreaterGreater':'\u2AA2\u0338','NotNestedLessLess':'\u2AA1\u0338','notni':'\u220C','notniva':'\u220C','notnivb':'\u22FE','notnivc':'\u22FD','NotPrecedes':'\u2280','NotPrecedesEqual':'\u2AAF\u0338','NotPrecedesSlantEqual':'\u22E0','NotReverseElement':'\u220C','NotRightTriangle':'\u22EB','NotRightTriangleBar':'\u29D0\u0338','NotRightTriangleEqual':'\u22ED','NotSquareSubset':'\u228F\u0338','NotSquareSubsetEqual':'\u22E2','NotSquareSuperset':'\u2290\u0338','NotSquareSupersetEqual':'\u22E3','NotSubset':'\u2282\u20D2','NotSubsetEqual':'\u2288','NotSucceeds':'\u2281','NotSucceedsEqual':'\u2AB0\u0338','NotSucceedsSlantEqual':'\u22E1','NotSucceedsTilde':'\u227F\u0338','NotSuperset':'\u2283\u20D2','NotSupersetEqual':'\u2289','NotTilde':'\u2241','NotTildeEqual':'\u2244','NotTildeFullEqual':'\u2247','NotTildeTilde':'\u2249','NotVerticalBar':'\u2224','npar':'\u2226','nparallel':'\u2226','nparsl':'\u2AFD\u20E5','npart':'\u2202\u0338','npolint':'\u2A14','npr':'\u2280','nprcue':'\u22E0','npre':'\u2AAF\u0338','nprec':'\u2280','npreceq':'\u2AAF\u0338','nrarr':'\u219B','nrArr':'\u21CF','nrarrc':'\u2933\u0338','nrarrw':'\u219D\u0338','nrightarrow':'\u219B','nRightarrow':'\u21CF','nrtri':'\u22EB','nrtrie':'\u22ED','nsc':'\u2281','nsccue':'\u22E1','nsce':'\u2AB0\u0338','nscr':'\uD835\uDCC3','Nscr':'\uD835\uDCA9','nshortmid':'\u2224','nshortparallel':'\u2226','nsim':'\u2241','nsime':'\u2244','nsimeq':'\u2244','nsmid':'\u2224','nspar':'\u2226','nsqsube':'\u22E2','nsqsupe':'\u22E3','nsub':'\u2284','nsube':'\u2288','nsubE':'\u2AC5\u0338','nsubset':'\u2282\u20D2','nsubseteq':'\u2288','nsubseteqq':'\u2AC5\u0338','nsucc':'\u2281','nsucceq':'\u2AB0\u0338','nsup':'\u2285','nsupe':'\u2289','nsupE':'\u2AC6\u0338','nsupset':'\u2283\u20D2','nsupseteq':'\u2289','nsupseteqq':'\u2AC6\u0338','ntgl':'\u2279','ntilde':'\xF1','Ntilde':'\xD1','ntlg':'\u2278','ntriangleleft':'\u22EA','ntrianglelefteq':'\u22EC','ntriangleright':'\u22EB','ntrianglerighteq':'\u22ED','nu':'\u03BD','Nu':'\u039D','num':'#','numero':'\u2116','numsp':'\u2007','nvap':'\u224D\u20D2','nvdash':'\u22AC','nvDash':'\u22AD','nVdash':'\u22AE','nVDash':'\u22AF','nvge':'\u2265\u20D2','nvgt':'>\u20D2','nvHarr':'\u2904','nvinfin':'\u29DE','nvlArr':'\u2902','nvle':'\u2264\u20D2','nvlt':'<\u20D2','nvltrie':'\u22B4\u20D2','nvrArr':'\u2903','nvrtrie':'\u22B5\u20D2','nvsim':'\u223C\u20D2','nwarhk':'\u2923','nwarr':'\u2196','nwArr':'\u21D6','nwarrow':'\u2196','nwnear':'\u2927','oacute':'\xF3','Oacute':'\xD3','oast':'\u229B','ocir':'\u229A','ocirc':'\xF4','Ocirc':'\xD4','ocy':'\u043E','Ocy':'\u041E','odash':'\u229D','odblac':'\u0151','Odblac':'\u0150','odiv':'\u2A38','odot':'\u2299','odsold':'\u29BC','oelig':'\u0153','OElig':'\u0152','ofcir':'\u29BF','ofr':'\uD835\uDD2C','Ofr':'\uD835\uDD12','ogon':'\u02DB','ograve':'\xF2','Ograve':'\xD2','ogt':'\u29C1','ohbar':'\u29B5','ohm':'\u03A9','oint':'\u222E','olarr':'\u21BA','olcir':'\u29BE','olcross':'\u29BB','oline':'\u203E','olt':'\u29C0','omacr':'\u014D','Omacr':'\u014C','omega':'\u03C9','Omega':'\u03A9','omicron':'\u03BF','Omicron':'\u039F','omid':'\u29B6','ominus':'\u2296','oopf':'\uD835\uDD60','Oopf':'\uD835\uDD46','opar':'\u29B7','OpenCurlyDoubleQuote':'\u201C','OpenCurlyQuote':'\u2018','operp':'\u29B9','oplus':'\u2295','or':'\u2228','Or':'\u2A54','orarr':'\u21BB','ord':'\u2A5D','order':'\u2134','orderof':'\u2134','ordf':'\xAA','ordm':'\xBA','origof':'\u22B6','oror':'\u2A56','orslope':'\u2A57','orv':'\u2A5B','oS':'\u24C8','oscr':'\u2134','Oscr':'\uD835\uDCAA','oslash':'\xF8','Oslash':'\xD8','osol':'\u2298','otilde':'\xF5','Otilde':'\xD5','otimes':'\u2297','Otimes':'\u2A37','otimesas':'\u2A36','ouml':'\xF6','Ouml':'\xD6','ovbar':'\u233D','OverBar':'\u203E','OverBrace':'\u23DE','OverBracket':'\u23B4','OverParenthesis':'\u23DC','par':'\u2225','para':'\xB6','parallel':'\u2225','parsim':'\u2AF3','parsl':'\u2AFD','part':'\u2202','PartialD':'\u2202','pcy':'\u043F','Pcy':'\u041F','percnt':'%','period':'.','permil':'\u2030','perp':'\u22A5','pertenk':'\u2031','pfr':'\uD835\uDD2D','Pfr':'\uD835\uDD13','phi':'\u03C6','Phi':'\u03A6','phiv':'\u03D5','phmmat':'\u2133','phone':'\u260E','pi':'\u03C0','Pi':'\u03A0','pitchfork':'\u22D4','piv':'\u03D6','planck':'\u210F','planckh':'\u210E','plankv':'\u210F','plus':'+','plusacir':'\u2A23','plusb':'\u229E','pluscir':'\u2A22','plusdo':'\u2214','plusdu':'\u2A25','pluse':'\u2A72','PlusMinus':'\xB1','plusmn':'\xB1','plussim':'\u2A26','plustwo':'\u2A27','pm':'\xB1','Poincareplane':'\u210C','pointint':'\u2A15','popf':'\uD835\uDD61','Popf':'\u2119','pound':'\xA3','pr':'\u227A','Pr':'\u2ABB','prap':'\u2AB7','prcue':'\u227C','pre':'\u2AAF','prE':'\u2AB3','prec':'\u227A','precapprox':'\u2AB7','preccurlyeq':'\u227C','Precedes':'\u227A','PrecedesEqual':'\u2AAF','PrecedesSlantEqual':'\u227C','PrecedesTilde':'\u227E','preceq':'\u2AAF','precnapprox':'\u2AB9','precneqq':'\u2AB5','precnsim':'\u22E8','precsim':'\u227E','prime':'\u2032','Prime':'\u2033','primes':'\u2119','prnap':'\u2AB9','prnE':'\u2AB5','prnsim':'\u22E8','prod':'\u220F','Product':'\u220F','profalar':'\u232E','profline':'\u2312','profsurf':'\u2313','prop':'\u221D','Proportion':'\u2237','Proportional':'\u221D','propto':'\u221D','prsim':'\u227E','prurel':'\u22B0','pscr':'\uD835\uDCC5','Pscr':'\uD835\uDCAB','psi':'\u03C8','Psi':'\u03A8','puncsp':'\u2008','qfr':'\uD835\uDD2E','Qfr':'\uD835\uDD14','qint':'\u2A0C','qopf':'\uD835\uDD62','Qopf':'\u211A','qprime':'\u2057','qscr':'\uD835\uDCC6','Qscr':'\uD835\uDCAC','quaternions':'\u210D','quatint':'\u2A16','quest':'?','questeq':'\u225F','quot':'"','QUOT':'"','rAarr':'\u21DB','race':'\u223D\u0331','racute':'\u0155','Racute':'\u0154','radic':'\u221A','raemptyv':'\u29B3','rang':'\u27E9','Rang':'\u27EB','rangd':'\u2992','range':'\u29A5','rangle':'\u27E9','raquo':'\xBB','rarr':'\u2192','rArr':'\u21D2','Rarr':'\u21A0','rarrap':'\u2975','rarrb':'\u21E5','rarrbfs':'\u2920','rarrc':'\u2933','rarrfs':'\u291E','rarrhk':'\u21AA','rarrlp':'\u21AC','rarrpl':'\u2945','rarrsim':'\u2974','rarrtl':'\u21A3','Rarrtl':'\u2916','rarrw':'\u219D','ratail':'\u291A','rAtail':'\u291C','ratio':'\u2236','rationals':'\u211A','rbarr':'\u290D','rBarr':'\u290F','RBarr':'\u2910','rbbrk':'\u2773','rbrace':'}','rbrack':']','rbrke':'\u298C','rbrksld':'\u298E','rbrkslu':'\u2990','rcaron':'\u0159','Rcaron':'\u0158','rcedil':'\u0157','Rcedil':'\u0156','rceil':'\u2309','rcub':'}','rcy':'\u0440','Rcy':'\u0420','rdca':'\u2937','rdldhar':'\u2969','rdquo':'\u201D','rdquor':'\u201D','rdsh':'\u21B3','Re':'\u211C','real':'\u211C','realine':'\u211B','realpart':'\u211C','reals':'\u211D','rect':'\u25AD','reg':'\xAE','REG':'\xAE','ReverseElement':'\u220B','ReverseEquilibrium':'\u21CB','ReverseUpEquilibrium':'\u296F','rfisht':'\u297D','rfloor':'\u230B','rfr':'\uD835\uDD2F','Rfr':'\u211C','rHar':'\u2964','rhard':'\u21C1','rharu':'\u21C0','rharul':'\u296C','rho':'\u03C1','Rho':'\u03A1','rhov':'\u03F1','RightAngleBracket':'\u27E9','rightarrow':'\u2192','Rightarrow':'\u21D2','RightArrow':'\u2192','RightArrowBar':'\u21E5','RightArrowLeftArrow':'\u21C4','rightarrowtail':'\u21A3','RightCeiling':'\u2309','RightDoubleBracket':'\u27E7','RightDownTeeVector':'\u295D','RightDownVector':'\u21C2','RightDownVectorBar':'\u2955','RightFloor':'\u230B','rightharpoondown':'\u21C1','rightharpoonup':'\u21C0','rightleftarrows':'\u21C4','rightleftharpoons':'\u21CC','rightrightarrows':'\u21C9','rightsquigarrow':'\u219D','RightTee':'\u22A2','RightTeeArrow':'\u21A6','RightTeeVector':'\u295B','rightthreetimes':'\u22CC','RightTriangle':'\u22B3','RightTriangleBar':'\u29D0','RightTriangleEqual':'\u22B5','RightUpDownVector':'\u294F','RightUpTeeVector':'\u295C','RightUpVector':'\u21BE','RightUpVectorBar':'\u2954','RightVector':'\u21C0','RightVectorBar':'\u2953','ring':'\u02DA','risingdotseq':'\u2253','rlarr':'\u21C4','rlhar':'\u21CC','rlm':'\u200F','rmoust':'\u23B1','rmoustache':'\u23B1','rnmid':'\u2AEE','roang':'\u27ED','roarr':'\u21FE','robrk':'\u27E7','ropar':'\u2986','ropf':'\uD835\uDD63','Ropf':'\u211D','roplus':'\u2A2E','rotimes':'\u2A35','RoundImplies':'\u2970','rpar':')','rpargt':'\u2994','rppolint':'\u2A12','rrarr':'\u21C9','Rrightarrow':'\u21DB','rsaquo':'\u203A','rscr':'\uD835\uDCC7','Rscr':'\u211B','rsh':'\u21B1','Rsh':'\u21B1','rsqb':']','rsquo':'\u2019','rsquor':'\u2019','rthree':'\u22CC','rtimes':'\u22CA','rtri':'\u25B9','rtrie':'\u22B5','rtrif':'\u25B8','rtriltri':'\u29CE','RuleDelayed':'\u29F4','ruluhar':'\u2968','rx':'\u211E','sacute':'\u015B','Sacute':'\u015A','sbquo':'\u201A','sc':'\u227B','Sc':'\u2ABC','scap':'\u2AB8','scaron':'\u0161','Scaron':'\u0160','sccue':'\u227D','sce':'\u2AB0','scE':'\u2AB4','scedil':'\u015F','Scedil':'\u015E','scirc':'\u015D','Scirc':'\u015C','scnap':'\u2ABA','scnE':'\u2AB6','scnsim':'\u22E9','scpolint':'\u2A13','scsim':'\u227F','scy':'\u0441','Scy':'\u0421','sdot':'\u22C5','sdotb':'\u22A1','sdote':'\u2A66','searhk':'\u2925','searr':'\u2198','seArr':'\u21D8','searrow':'\u2198','sect':'\xA7','semi':';','seswar':'\u2929','setminus':'\u2216','setmn':'\u2216','sext':'\u2736','sfr':'\uD835\uDD30','Sfr':'\uD835\uDD16','sfrown':'\u2322','sharp':'\u266F','shchcy':'\u0449','SHCHcy':'\u0429','shcy':'\u0448','SHcy':'\u0428','ShortDownArrow':'\u2193','ShortLeftArrow':'\u2190','shortmid':'\u2223','shortparallel':'\u2225','ShortRightArrow':'\u2192','ShortUpArrow':'\u2191','shy':'\xAD','sigma':'\u03C3','Sigma':'\u03A3','sigmaf':'\u03C2','sigmav':'\u03C2','sim':'\u223C','simdot':'\u2A6A','sime':'\u2243','simeq':'\u2243','simg':'\u2A9E','simgE':'\u2AA0','siml':'\u2A9D','simlE':'\u2A9F','simne':'\u2246','simplus':'\u2A24','simrarr':'\u2972','slarr':'\u2190','SmallCircle':'\u2218','smallsetminus':'\u2216','smashp':'\u2A33','smeparsl':'\u29E4','smid':'\u2223','smile':'\u2323','smt':'\u2AAA','smte':'\u2AAC','smtes':'\u2AAC\uFE00','softcy':'\u044C','SOFTcy':'\u042C','sol':'/','solb':'\u29C4','solbar':'\u233F','sopf':'\uD835\uDD64','Sopf':'\uD835\uDD4A','spades':'\u2660','spadesuit':'\u2660','spar':'\u2225','sqcap':'\u2293','sqcaps':'\u2293\uFE00','sqcup':'\u2294','sqcups':'\u2294\uFE00','Sqrt':'\u221A','sqsub':'\u228F','sqsube':'\u2291','sqsubset':'\u228F','sqsubseteq':'\u2291','sqsup':'\u2290','sqsupe':'\u2292','sqsupset':'\u2290','sqsupseteq':'\u2292','squ':'\u25A1','square':'\u25A1','Square':'\u25A1','SquareIntersection':'\u2293','SquareSubset':'\u228F','SquareSubsetEqual':'\u2291','SquareSuperset':'\u2290','SquareSupersetEqual':'\u2292','SquareUnion':'\u2294','squarf':'\u25AA','squf':'\u25AA','srarr':'\u2192','sscr':'\uD835\uDCC8','Sscr':'\uD835\uDCAE','ssetmn':'\u2216','ssmile':'\u2323','sstarf':'\u22C6','star':'\u2606','Star':'\u22C6','starf':'\u2605','straightepsilon':'\u03F5','straightphi':'\u03D5','strns':'\xAF','sub':'\u2282','Sub':'\u22D0','subdot':'\u2ABD','sube':'\u2286','subE':'\u2AC5','subedot':'\u2AC3','submult':'\u2AC1','subne':'\u228A','subnE':'\u2ACB','subplus':'\u2ABF','subrarr':'\u2979','subset':'\u2282','Subset':'\u22D0','subseteq':'\u2286','subseteqq':'\u2AC5','SubsetEqual':'\u2286','subsetneq':'\u228A','subsetneqq':'\u2ACB','subsim':'\u2AC7','subsub':'\u2AD5','subsup':'\u2AD3','succ':'\u227B','succapprox':'\u2AB8','succcurlyeq':'\u227D','Succeeds':'\u227B','SucceedsEqual':'\u2AB0','SucceedsSlantEqual':'\u227D','SucceedsTilde':'\u227F','succeq':'\u2AB0','succnapprox':'\u2ABA','succneqq':'\u2AB6','succnsim':'\u22E9','succsim':'\u227F','SuchThat':'\u220B','sum':'\u2211','Sum':'\u2211','sung':'\u266A','sup':'\u2283','Sup':'\u22D1','sup1':'\xB9','sup2':'\xB2','sup3':'\xB3','supdot':'\u2ABE','supdsub':'\u2AD8','supe':'\u2287','supE':'\u2AC6','supedot':'\u2AC4','Superset':'\u2283','SupersetEqual':'\u2287','suphsol':'\u27C9','suphsub':'\u2AD7','suplarr':'\u297B','supmult':'\u2AC2','supne':'\u228B','supnE':'\u2ACC','supplus':'\u2AC0','supset':'\u2283','Supset':'\u22D1','supseteq':'\u2287','supseteqq':'\u2AC6','supsetneq':'\u228B','supsetneqq':'\u2ACC','supsim':'\u2AC8','supsub':'\u2AD4','supsup':'\u2AD6','swarhk':'\u2926','swarr':'\u2199','swArr':'\u21D9','swarrow':'\u2199','swnwar':'\u292A','szlig':'\xDF','Tab':'\t','target':'\u2316','tau':'\u03C4','Tau':'\u03A4','tbrk':'\u23B4','tcaron':'\u0165','Tcaron':'\u0164','tcedil':'\u0163','Tcedil':'\u0162','tcy':'\u0442','Tcy':'\u0422','tdot':'\u20DB','telrec':'\u2315','tfr':'\uD835\uDD31','Tfr':'\uD835\uDD17','there4':'\u2234','therefore':'\u2234','Therefore':'\u2234','theta':'\u03B8','Theta':'\u0398','thetasym':'\u03D1','thetav':'\u03D1','thickapprox':'\u2248','thicksim':'\u223C','ThickSpace':'\u205F\u200A','thinsp':'\u2009','ThinSpace':'\u2009','thkap':'\u2248','thksim':'\u223C','thorn':'\xFE','THORN':'\xDE','tilde':'\u02DC','Tilde':'\u223C','TildeEqual':'\u2243','TildeFullEqual':'\u2245','TildeTilde':'\u2248','times':'\xD7','timesb':'\u22A0','timesbar':'\u2A31','timesd':'\u2A30','tint':'\u222D','toea':'\u2928','top':'\u22A4','topbot':'\u2336','topcir':'\u2AF1','topf':'\uD835\uDD65','Topf':'\uD835\uDD4B','topfork':'\u2ADA','tosa':'\u2929','tprime':'\u2034','trade':'\u2122','TRADE':'\u2122','triangle':'\u25B5','triangledown':'\u25BF','triangleleft':'\u25C3','trianglelefteq':'\u22B4','triangleq':'\u225C','triangleright':'\u25B9','trianglerighteq':'\u22B5','tridot':'\u25EC','trie':'\u225C','triminus':'\u2A3A','TripleDot':'\u20DB','triplus':'\u2A39','trisb':'\u29CD','tritime':'\u2A3B','trpezium':'\u23E2','tscr':'\uD835\uDCC9','Tscr':'\uD835\uDCAF','tscy':'\u0446','TScy':'\u0426','tshcy':'\u045B','TSHcy':'\u040B','tstrok':'\u0167','Tstrok':'\u0166','twixt':'\u226C','twoheadleftarrow':'\u219E','twoheadrightarrow':'\u21A0','uacute':'\xFA','Uacute':'\xDA','uarr':'\u2191','uArr':'\u21D1','Uarr':'\u219F','Uarrocir':'\u2949','ubrcy':'\u045E','Ubrcy':'\u040E','ubreve':'\u016D','Ubreve':'\u016C','ucirc':'\xFB','Ucirc':'\xDB','ucy':'\u0443','Ucy':'\u0423','udarr':'\u21C5','udblac':'\u0171','Udblac':'\u0170','udhar':'\u296E','ufisht':'\u297E','ufr':'\uD835\uDD32','Ufr':'\uD835\uDD18','ugrave':'\xF9','Ugrave':'\xD9','uHar':'\u2963','uharl':'\u21BF','uharr':'\u21BE','uhblk':'\u2580','ulcorn':'\u231C','ulcorner':'\u231C','ulcrop':'\u230F','ultri':'\u25F8','umacr':'\u016B','Umacr':'\u016A','uml':'\xA8','UnderBar':'_','UnderBrace':'\u23DF','UnderBracket':'\u23B5','UnderParenthesis':'\u23DD','Union':'\u22C3','UnionPlus':'\u228E','uogon':'\u0173','Uogon':'\u0172','uopf':'\uD835\uDD66','Uopf':'\uD835\uDD4C','uparrow':'\u2191','Uparrow':'\u21D1','UpArrow':'\u2191','UpArrowBar':'\u2912','UpArrowDownArrow':'\u21C5','updownarrow':'\u2195','Updownarrow':'\u21D5','UpDownArrow':'\u2195','UpEquilibrium':'\u296E','upharpoonleft':'\u21BF','upharpoonright':'\u21BE','uplus':'\u228E','UpperLeftArrow':'\u2196','UpperRightArrow':'\u2197','upsi':'\u03C5','Upsi':'\u03D2','upsih':'\u03D2','upsilon':'\u03C5','Upsilon':'\u03A5','UpTee':'\u22A5','UpTeeArrow':'\u21A5','upuparrows':'\u21C8','urcorn':'\u231D','urcorner':'\u231D','urcrop':'\u230E','uring':'\u016F','Uring':'\u016E','urtri':'\u25F9','uscr':'\uD835\uDCCA','Uscr':'\uD835\uDCB0','utdot':'\u22F0','utilde':'\u0169','Utilde':'\u0168','utri':'\u25B5','utrif':'\u25B4','uuarr':'\u21C8','uuml':'\xFC','Uuml':'\xDC','uwangle':'\u29A7','vangrt':'\u299C','varepsilon':'\u03F5','varkappa':'\u03F0','varnothing':'\u2205','varphi':'\u03D5','varpi':'\u03D6','varpropto':'\u221D','varr':'\u2195','vArr':'\u21D5','varrho':'\u03F1','varsigma':'\u03C2','varsubsetneq':'\u228A\uFE00','varsubsetneqq':'\u2ACB\uFE00','varsupsetneq':'\u228B\uFE00','varsupsetneqq':'\u2ACC\uFE00','vartheta':'\u03D1','vartriangleleft':'\u22B2','vartriangleright':'\u22B3','vBar':'\u2AE8','Vbar':'\u2AEB','vBarv':'\u2AE9','vcy':'\u0432','Vcy':'\u0412','vdash':'\u22A2','vDash':'\u22A8','Vdash':'\u22A9','VDash':'\u22AB','Vdashl':'\u2AE6','vee':'\u2228','Vee':'\u22C1','veebar':'\u22BB','veeeq':'\u225A','vellip':'\u22EE','verbar':'|','Verbar':'\u2016','vert':'|','Vert':'\u2016','VerticalBar':'\u2223','VerticalLine':'|','VerticalSeparator':'\u2758','VerticalTilde':'\u2240','VeryThinSpace':'\u200A','vfr':'\uD835\uDD33','Vfr':'\uD835\uDD19','vltri':'\u22B2','vnsub':'\u2282\u20D2','vnsup':'\u2283\u20D2','vopf':'\uD835\uDD67','Vopf':'\uD835\uDD4D','vprop':'\u221D','vrtri':'\u22B3','vscr':'\uD835\uDCCB','Vscr':'\uD835\uDCB1','vsubne':'\u228A\uFE00','vsubnE':'\u2ACB\uFE00','vsupne':'\u228B\uFE00','vsupnE':'\u2ACC\uFE00','Vvdash':'\u22AA','vzigzag':'\u299A','wcirc':'\u0175','Wcirc':'\u0174','wedbar':'\u2A5F','wedge':'\u2227','Wedge':'\u22C0','wedgeq':'\u2259','weierp':'\u2118','wfr':'\uD835\uDD34','Wfr':'\uD835\uDD1A','wopf':'\uD835\uDD68','Wopf':'\uD835\uDD4E','wp':'\u2118','wr':'\u2240','wreath':'\u2240','wscr':'\uD835\uDCCC','Wscr':'\uD835\uDCB2','xcap':'\u22C2','xcirc':'\u25EF','xcup':'\u22C3','xdtri':'\u25BD','xfr':'\uD835\uDD35','Xfr':'\uD835\uDD1B','xharr':'\u27F7','xhArr':'\u27FA','xi':'\u03BE','Xi':'\u039E','xlarr':'\u27F5','xlArr':'\u27F8','xmap':'\u27FC','xnis':'\u22FB','xodot':'\u2A00','xopf':'\uD835\uDD69','Xopf':'\uD835\uDD4F','xoplus':'\u2A01','xotime':'\u2A02','xrarr':'\u27F6','xrArr':'\u27F9','xscr':'\uD835\uDCCD','Xscr':'\uD835\uDCB3','xsqcup':'\u2A06','xuplus':'\u2A04','xutri':'\u25B3','xvee':'\u22C1','xwedge':'\u22C0','yacute':'\xFD','Yacute':'\xDD','yacy':'\u044F','YAcy':'\u042F','ycirc':'\u0177','Ycirc':'\u0176','ycy':'\u044B','Ycy':'\u042B','yen':'\xA5','yfr':'\uD835\uDD36','Yfr':'\uD835\uDD1C','yicy':'\u0457','YIcy':'\u0407','yopf':'\uD835\uDD6A','Yopf':'\uD835\uDD50','yscr':'\uD835\uDCCE','Yscr':'\uD835\uDCB4','yucy':'\u044E','YUcy':'\u042E','yuml':'\xFF','Yuml':'\u0178','zacute':'\u017A','Zacute':'\u0179','zcaron':'\u017E','Zcaron':'\u017D','zcy':'\u0437','Zcy':'\u0417','zdot':'\u017C','Zdot':'\u017B','zeetrf':'\u2128','ZeroWidthSpace':'\u200B','zeta':'\u03B6','Zeta':'\u0396','zfr':'\uD835\uDD37','Zfr':'\u2128','zhcy':'\u0436','ZHcy':'\u0416','zigrarr':'\u21DD','zopf':'\uD835\uDD6B','Zopf':'\u2124','zscr':'\uD835\uDCCF','Zscr':'\uD835\uDCB5','zwj':'\u200D','zwnj':'\u200C'};
		var decodeMapLegacy = {'aacute':'\xE1','Aacute':'\xC1','acirc':'\xE2','Acirc':'\xC2','acute':'\xB4','aelig':'\xE6','AElig':'\xC6','agrave':'\xE0','Agrave':'\xC0','amp':'&','AMP':'&','aring':'\xE5','Aring':'\xC5','atilde':'\xE3','Atilde':'\xC3','auml':'\xE4','Auml':'\xC4','brvbar':'\xA6','ccedil':'\xE7','Ccedil':'\xC7','cedil':'\xB8','cent':'\xA2','copy':'\xA9','COPY':'\xA9','curren':'\xA4','deg':'\xB0','divide':'\xF7','eacute':'\xE9','Eacute':'\xC9','ecirc':'\xEA','Ecirc':'\xCA','egrave':'\xE8','Egrave':'\xC8','eth':'\xF0','ETH':'\xD0','euml':'\xEB','Euml':'\xCB','frac12':'\xBD','frac14':'\xBC','frac34':'\xBE','gt':'>','GT':'>','iacute':'\xED','Iacute':'\xCD','icirc':'\xEE','Icirc':'\xCE','iexcl':'\xA1','igrave':'\xEC','Igrave':'\xCC','iquest':'\xBF','iuml':'\xEF','Iuml':'\xCF','laquo':'\xAB','lt':'<','LT':'<','macr':'\xAF','micro':'\xB5','middot':'\xB7','nbsp':'\xA0','not':'\xAC','ntilde':'\xF1','Ntilde':'\xD1','oacute':'\xF3','Oacute':'\xD3','ocirc':'\xF4','Ocirc':'\xD4','ograve':'\xF2','Ograve':'\xD2','ordf':'\xAA','ordm':'\xBA','oslash':'\xF8','Oslash':'\xD8','otilde':'\xF5','Otilde':'\xD5','ouml':'\xF6','Ouml':'\xD6','para':'\xB6','plusmn':'\xB1','pound':'\xA3','quot':'"','QUOT':'"','raquo':'\xBB','reg':'\xAE','REG':'\xAE','sect':'\xA7','shy':'\xAD','sup1':'\xB9','sup2':'\xB2','sup3':'\xB3','szlig':'\xDF','thorn':'\xFE','THORN':'\xDE','times':'\xD7','uacute':'\xFA','Uacute':'\xDA','ucirc':'\xFB','Ucirc':'\xDB','ugrave':'\xF9','Ugrave':'\xD9','uml':'\xA8','uuml':'\xFC','Uuml':'\xDC','yacute':'\xFD','Yacute':'\xDD','yen':'\xA5','yuml':'\xFF'};
		var decodeMapNumeric = {'0':'\uFFFD','128':'\u20AC','130':'\u201A','131':'\u0192','132':'\u201E','133':'\u2026','134':'\u2020','135':'\u2021','136':'\u02C6','137':'\u2030','138':'\u0160','139':'\u2039','140':'\u0152','142':'\u017D','145':'\u2018','146':'\u2019','147':'\u201C','148':'\u201D','149':'\u2022','150':'\u2013','151':'\u2014','152':'\u02DC','153':'\u2122','154':'\u0161','155':'\u203A','156':'\u0153','158':'\u017E','159':'\u0178'};
		var invalidReferenceCodePoints = [1,2,3,4,5,6,7,8,11,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,64976,64977,64978,64979,64980,64981,64982,64983,64984,64985,64986,64987,64988,64989,64990,64991,64992,64993,64994,64995,64996,64997,64998,64999,65000,65001,65002,65003,65004,65005,65006,65007,65534,65535,131070,131071,196606,196607,262142,262143,327678,327679,393214,393215,458750,458751,524286,524287,589822,589823,655358,655359,720894,720895,786430,786431,851966,851967,917502,917503,983038,983039,1048574,1048575,1114110,1114111];

		/*--------------------------------------------------------------------------*/

		var stringFromCharCode = String.fromCharCode;

		var object = {};
		var hasOwnProperty = object.hasOwnProperty;
		var has = function(object, propertyName) {
			return hasOwnProperty.call(object, propertyName);
		};

		var contains = function(array, value) {
			var index = -1;
			var length = array.length;
			while (++index < length) {
				if (array[index] == value) {
					return true;
				}
			}
			return false;
		};

		var merge = function(options, defaults) {
			if (!options) {
				return defaults;
			}
			var result = {};
			var key;
			for (key in defaults) {
				// A `hasOwnProperty` check is not needed here, since only recognized
				// option names are used anyway. Any others are ignored.
				result[key] = has(options, key) ? options[key] : defaults[key];
			}
			return result;
		};

		// Modified version of `ucs2encode`; see https://mths.be/punycode.
		var codePointToSymbol = function(codePoint, strict) {
			var output = '';
			if ((codePoint >= 0xD800 && codePoint <= 0xDFFF) || codePoint > 0x10FFFF) {
				// See issue #4:
				// “Otherwise, if the number is in the range 0xD800 to 0xDFFF or is
				// greater than 0x10FFFF, then this is a parse error. Return a U+FFFD
				// REPLACEMENT CHARACTER.”
				if (strict) {
					parseError('character reference outside the permissible Unicode range');
				}
				return '\uFFFD';
			}
			if (has(decodeMapNumeric, codePoint)) {
				if (strict) {
					parseError('disallowed character reference');
				}
				return decodeMapNumeric[codePoint];
			}
			if (strict && contains(invalidReferenceCodePoints, codePoint)) {
				parseError('disallowed character reference');
			}
			if (codePoint > 0xFFFF) {
				codePoint -= 0x10000;
				output += stringFromCharCode(codePoint >>> 10 & 0x3FF | 0xD800);
				codePoint = 0xDC00 | codePoint & 0x3FF;
			}
			output += stringFromCharCode(codePoint);
			return output;
		};

		var hexEscape = function(codePoint) {
			return '&#x' + codePoint.toString(16).toUpperCase() + ';';
		};

		var decEscape = function(codePoint) {
			return '&#' + codePoint + ';';
		};

		var parseError = function(message) {
			throw Error('Parse error: ' + message);
		};

		/*--------------------------------------------------------------------------*/

		var encode = function(string, options) {
			options = merge(options, encode.options);
			var strict = options.strict;
			if (strict && regexInvalidRawCodePoint.test(string)) {
				parseError('forbidden code point');
			}
			var encodeEverything = options.encodeEverything;
			var useNamedReferences = options.useNamedReferences;
			var allowUnsafeSymbols = options.allowUnsafeSymbols;
			var escapeCodePoint = options.decimal ? decEscape : hexEscape;

			var escapeBmpSymbol = function(symbol) {
				return escapeCodePoint(symbol.charCodeAt(0));
			};

			if (encodeEverything) {
				// Encode ASCII symbols.
				string = string.replace(regexAsciiWhitelist, function(symbol) {
					// Use named references if requested & possible.
					if (useNamedReferences && has(encodeMap, symbol)) {
						return '&' + encodeMap[symbol] + ';';
					}
					return escapeBmpSymbol(symbol);
				});
				// Shorten a few escapes that represent two symbols, of which at least one
				// is within the ASCII range.
				if (useNamedReferences) {
					string = string
						.replace(/&gt;\u20D2/g, '&nvgt;')
						.replace(/&lt;\u20D2/g, '&nvlt;')
						.replace(/&#x66;&#x6A;/g, '&fjlig;');
				}
				// Encode non-ASCII symbols.
				if (useNamedReferences) {
					// Encode non-ASCII symbols that can be replaced with a named reference.
					string = string.replace(regexEncodeNonAscii, function(string) {
						// Note: there is no need to check `has(encodeMap, string)` here.
						return '&' + encodeMap[string] + ';';
					});
				}
				// Note: any remaining non-ASCII symbols are handled outside of the `if`.
			} else if (useNamedReferences) {
				// Apply named character references.
				// Encode `<>"'&` using named character references.
				if (!allowUnsafeSymbols) {
					string = string.replace(regexEscape, function(string) {
						return '&' + encodeMap[string] + ';'; // no need to check `has()` here
					});
				}
				// Shorten escapes that represent two symbols, of which at least one is
				// `<>"'&`.
				string = string
					.replace(/&gt;\u20D2/g, '&nvgt;')
					.replace(/&lt;\u20D2/g, '&nvlt;');
				// Encode non-ASCII symbols that can be replaced with a named reference.
				string = string.replace(regexEncodeNonAscii, function(string) {
					// Note: there is no need to check `has(encodeMap, string)` here.
					return '&' + encodeMap[string] + ';';
				});
			} else if (!allowUnsafeSymbols) {
				// Encode `<>"'&` using hexadecimal escapes, now that they’re not handled
				// using named character references.
				string = string.replace(regexEscape, escapeBmpSymbol);
			}
			return string
				// Encode astral symbols.
				.replace(regexAstralSymbols, function($0) {
					// https://mathiasbynens.be/notes/javascript-encoding#surrogate-formulae
					var high = $0.charCodeAt(0);
					var low = $0.charCodeAt(1);
					var codePoint = (high - 0xD800) * 0x400 + low - 0xDC00 + 0x10000;
					return escapeCodePoint(codePoint);
				})
				// Encode any remaining BMP symbols that are not printable ASCII symbols
				// using a hexadecimal escape.
				.replace(regexBmpWhitelist, escapeBmpSymbol);
		};
		// Expose default options (so they can be overridden globally).
		encode.options = {
			'allowUnsafeSymbols': false,
			'encodeEverything': false,
			'strict': false,
			'useNamedReferences': false,
			'decimal' : false
		};

		var decode = function(html, options) {
			options = merge(options, decode.options);
			var strict = options.strict;
			if (strict && regexInvalidEntity.test(html)) {
				parseError('malformed character reference');
			}
			return html.replace(regexDecode, function($0, $1, $2, $3, $4, $5, $6, $7, $8) {
				var codePoint;
				var semicolon;
				var decDigits;
				var hexDigits;
				var reference;
				var next;

				if ($1) {
					reference = $1;
					// Note: there is no need to check `has(decodeMap, reference)`.
					return decodeMap[reference];
				}

				if ($2) {
					// Decode named character references without trailing `;`, e.g. `&amp`.
					// This is only a parse error if it gets converted to `&`, or if it is
					// followed by `=` in an attribute context.
					reference = $2;
					next = $3;
					if (next && options.isAttributeValue) {
						if (strict && next == '=') {
							parseError('`&` did not start a character reference');
						}
						return $0;
					} else {
						if (strict) {
							parseError(
								'named character reference was not terminated by a semicolon'
							);
						}
						// Note: there is no need to check `has(decodeMapLegacy, reference)`.
						return decodeMapLegacy[reference] + (next || '');
					}
				}

				if ($4) {
					// Decode decimal escapes, e.g. `&#119558;`.
					decDigits = $4;
					semicolon = $5;
					if (strict && !semicolon) {
						parseError('character reference was not terminated by a semicolon');
					}
					codePoint = parseInt(decDigits, 10);
					return codePointToSymbol(codePoint, strict);
				}

				if ($6) {
					// Decode hexadecimal escapes, e.g. `&#x1D306;`.
					hexDigits = $6;
					semicolon = $7;
					if (strict && !semicolon) {
						parseError('character reference was not terminated by a semicolon');
					}
					codePoint = parseInt(hexDigits, 16);
					return codePointToSymbol(codePoint, strict);
				}

				// If we’re still here, `if ($7)` is implied; it’s an ambiguous
				// ampersand for sure. https://mths.be/notes/ambiguous-ampersands
				if (strict) {
					parseError(
						'named character reference was not terminated by a semicolon'
					);
				}
				return $0;
			});
		};
		// Expose default options (so they can be overridden globally).
		decode.options = {
			'isAttributeValue': false,
			'strict': false
		};

		var escape = function(string) {
			return string.replace(regexEscape, function($0) {
				// Note: there is no need to check `has(escapeMap, $0)` here.
				return escapeMap[$0];
			});
		};

		/*--------------------------------------------------------------------------*/

		var he = {
			'version': '1.2.0',
			'encode': encode,
			'decode': decode,
			'escape': escape,
			'unescape': decode
		};

		// Some AMD build optimizers, like r.js, check for specific condition patterns
		// like the following:
		if (freeExports && !freeExports.nodeType) {
			if (freeModule) { // in Node.js, io.js, or RingoJS v0.8.0+
				freeModule.exports = he;
			} else { // in Narwhal or RingoJS v0.7.0-
				for (var key in he) {
					has(he, key) && (freeExports[key] = he[key]);
				}
			}
		} else { // in Rhino or a web browser
			root.he = he;
		}

	}(commonjsGlobal)); 
} (he$2, he$2.exports));

var heExports = he$2.exports;

var lib$3 = {};

var Parser$1 = {};

var Tokenizer$1 = {};

var decode_codepoint$1 = {};

var require$$0$2 = {
	"0": 65533,
	"128": 8364,
	"130": 8218,
	"131": 402,
	"132": 8222,
	"133": 8230,
	"134": 8224,
	"135": 8225,
	"136": 710,
	"137": 8240,
	"138": 352,
	"139": 8249,
	"140": 338,
	"142": 381,
	"145": 8216,
	"146": 8217,
	"147": 8220,
	"148": 8221,
	"149": 8226,
	"150": 8211,
	"151": 8212,
	"152": 732,
	"153": 8482,
	"154": 353,
	"155": 8250,
	"156": 339,
	"158": 382,
	"159": 376
};

var __importDefault$7 = (commonjsGlobal && commonjsGlobal.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(decode_codepoint$1, "__esModule", { value: true });
var decode_json_1$1 = __importDefault$7(require$$0$2);
// Adapted from https://github.com/mathiasbynens/he/blob/master/src/he.js#L94-L119
var fromCodePoint$1 = 
// eslint-disable-next-line @typescript-eslint/no-unnecessary-condition
String.fromCodePoint ||
    function (codePoint) {
        var output = "";
        if (codePoint > 0xffff) {
            codePoint -= 0x10000;
            output += String.fromCharCode(((codePoint >>> 10) & 0x3ff) | 0xd800);
            codePoint = 0xdc00 | (codePoint & 0x3ff);
        }
        output += String.fromCharCode(codePoint);
        return output;
    };
function decodeCodePoint$1(codePoint) {
    if ((codePoint >= 0xd800 && codePoint <= 0xdfff) || codePoint > 0x10ffff) {
        return "\uFFFD";
    }
    if (codePoint in decode_json_1$1.default) {
        codePoint = decode_json_1$1.default[codePoint];
    }
    return fromCodePoint$1(codePoint);
}
decode_codepoint$1.default = decodeCodePoint$1;

var Aacute$3 = "Á";
var aacute$3 = "á";
var Abreve$1 = "Ă";
var abreve$1 = "ă";
var ac$1 = "∾";
var acd$1 = "∿";
var acE$1 = "∾̳";
var Acirc$3 = "Â";
var acirc$3 = "â";
var acute$3 = "´";
var Acy$1 = "А";
var acy$1 = "а";
var AElig$3 = "Æ";
var aelig$3 = "æ";
var af$1 = "⁡";
var Afr$1 = "𝔄";
var afr$1 = "𝔞";
var Agrave$3 = "À";
var agrave$3 = "à";
var alefsym$1 = "ℵ";
var aleph$1 = "ℵ";
var Alpha$1 = "Α";
var alpha$1 = "α";
var Amacr$1 = "Ā";
var amacr$1 = "ā";
var amalg$1 = "⨿";
var amp$5 = "&";
var AMP$3 = "&";
var andand$1 = "⩕";
var And$1 = "⩓";
var and$1 = "∧";
var andd$1 = "⩜";
var andslope$1 = "⩘";
var andv$1 = "⩚";
var ang$1 = "∠";
var ange$1 = "⦤";
var angle$1 = "∠";
var angmsdaa$1 = "⦨";
var angmsdab$1 = "⦩";
var angmsdac$1 = "⦪";
var angmsdad$1 = "⦫";
var angmsdae$1 = "⦬";
var angmsdaf$1 = "⦭";
var angmsdag$1 = "⦮";
var angmsdah$1 = "⦯";
var angmsd$1 = "∡";
var angrt$1 = "∟";
var angrtvb$1 = "⊾";
var angrtvbd$1 = "⦝";
var angsph$1 = "∢";
var angst$1 = "Å";
var angzarr$1 = "⍼";
var Aogon$1 = "Ą";
var aogon$1 = "ą";
var Aopf$1 = "𝔸";
var aopf$1 = "𝕒";
var apacir$1 = "⩯";
var ap$1 = "≈";
var apE$1 = "⩰";
var ape$1 = "≊";
var apid$1 = "≋";
var apos$3 = "'";
var ApplyFunction$1 = "⁡";
var approx$1 = "≈";
var approxeq$1 = "≊";
var Aring$3 = "Å";
var aring$3 = "å";
var Ascr$1 = "𝒜";
var ascr$1 = "𝒶";
var Assign$1 = "≔";
var ast$1 = "*";
var asymp$1 = "≈";
var asympeq$1 = "≍";
var Atilde$3 = "Ã";
var atilde$3 = "ã";
var Auml$3 = "Ä";
var auml$3 = "ä";
var awconint$1 = "∳";
var awint$1 = "⨑";
var backcong$1 = "≌";
var backepsilon$1 = "϶";
var backprime$1 = "‵";
var backsim$1 = "∽";
var backsimeq$1 = "⋍";
var Backslash$1 = "∖";
var Barv$1 = "⫧";
var barvee$1 = "⊽";
var barwed$1 = "⌅";
var Barwed$1 = "⌆";
var barwedge$1 = "⌅";
var bbrk$1 = "⎵";
var bbrktbrk$1 = "⎶";
var bcong$1 = "≌";
var Bcy$1 = "Б";
var bcy$1 = "б";
var bdquo$1 = "„";
var becaus$1 = "∵";
var because$1 = "∵";
var Because$1 = "∵";
var bemptyv$1 = "⦰";
var bepsi$1 = "϶";
var bernou$1 = "ℬ";
var Bernoullis$1 = "ℬ";
var Beta$1 = "Β";
var beta$1 = "β";
var beth$1 = "ℶ";
var between$1 = "≬";
var Bfr$1 = "𝔅";
var bfr$1 = "𝔟";
var bigcap$1 = "⋂";
var bigcirc$1 = "◯";
var bigcup$1 = "⋃";
var bigodot$1 = "⨀";
var bigoplus$1 = "⨁";
var bigotimes$1 = "⨂";
var bigsqcup$1 = "⨆";
var bigstar$1 = "★";
var bigtriangledown$1 = "▽";
var bigtriangleup$1 = "△";
var biguplus$1 = "⨄";
var bigvee$1 = "⋁";
var bigwedge$1 = "⋀";
var bkarow$1 = "⤍";
var blacklozenge$1 = "⧫";
var blacksquare$1 = "▪";
var blacktriangle$1 = "▴";
var blacktriangledown$1 = "▾";
var blacktriangleleft$1 = "◂";
var blacktriangleright$1 = "▸";
var blank$1 = "␣";
var blk12$1 = "▒";
var blk14$1 = "░";
var blk34$1 = "▓";
var block$1 = "█";
var bne$1 = "=⃥";
var bnequiv$1 = "≡⃥";
var bNot$1 = "⫭";
var bnot$1 = "⌐";
var Bopf$1 = "𝔹";
var bopf$1 = "𝕓";
var bot$1 = "⊥";
var bottom$1 = "⊥";
var bowtie$1 = "⋈";
var boxbox$1 = "⧉";
var boxdl$1 = "┐";
var boxdL$1 = "╕";
var boxDl$1 = "╖";
var boxDL$1 = "╗";
var boxdr$1 = "┌";
var boxdR$1 = "╒";
var boxDr$1 = "╓";
var boxDR$1 = "╔";
var boxh$1 = "─";
var boxH$1 = "═";
var boxhd$1 = "┬";
var boxHd$1 = "╤";
var boxhD$1 = "╥";
var boxHD$1 = "╦";
var boxhu$1 = "┴";
var boxHu$1 = "╧";
var boxhU$1 = "╨";
var boxHU$1 = "╩";
var boxminus$1 = "⊟";
var boxplus$1 = "⊞";
var boxtimes$1 = "⊠";
var boxul$1 = "┘";
var boxuL$1 = "╛";
var boxUl$1 = "╜";
var boxUL$1 = "╝";
var boxur$1 = "└";
var boxuR$1 = "╘";
var boxUr$1 = "╙";
var boxUR$1 = "╚";
var boxv$1 = "│";
var boxV$1 = "║";
var boxvh$1 = "┼";
var boxvH$1 = "╪";
var boxVh$1 = "╫";
var boxVH$1 = "╬";
var boxvl$1 = "┤";
var boxvL$1 = "╡";
var boxVl$1 = "╢";
var boxVL$1 = "╣";
var boxvr$1 = "├";
var boxvR$1 = "╞";
var boxVr$1 = "╟";
var boxVR$1 = "╠";
var bprime$1 = "‵";
var breve$1 = "˘";
var Breve$1 = "˘";
var brvbar$3 = "¦";
var bscr$1 = "𝒷";
var Bscr$1 = "ℬ";
var bsemi$1 = "⁏";
var bsim$1 = "∽";
var bsime$1 = "⋍";
var bsolb$1 = "⧅";
var bsol$1 = "\\";
var bsolhsub$1 = "⟈";
var bull$1 = "•";
var bullet$1 = "•";
var bump$1 = "≎";
var bumpE$1 = "⪮";
var bumpe$1 = "≏";
var Bumpeq$1 = "≎";
var bumpeq$1 = "≏";
var Cacute$1 = "Ć";
var cacute$1 = "ć";
var capand$1 = "⩄";
var capbrcup$1 = "⩉";
var capcap$1 = "⩋";
var cap$1 = "∩";
var Cap$1 = "⋒";
var capcup$1 = "⩇";
var capdot$1 = "⩀";
var CapitalDifferentialD$1 = "ⅅ";
var caps$1 = "∩︀";
var caret$2 = "⁁";
var caron$1 = "ˇ";
var Cayleys$1 = "ℭ";
var ccaps$1 = "⩍";
var Ccaron$1 = "Č";
var ccaron$1 = "č";
var Ccedil$3 = "Ç";
var ccedil$3 = "ç";
var Ccirc$1 = "Ĉ";
var ccirc$1 = "ĉ";
var Cconint$1 = "∰";
var ccups$1 = "⩌";
var ccupssm$1 = "⩐";
var Cdot$1 = "Ċ";
var cdot$1 = "ċ";
var cedil$3 = "¸";
var Cedilla$1 = "¸";
var cemptyv$1 = "⦲";
var cent$3 = "¢";
var centerdot$1 = "·";
var CenterDot$1 = "·";
var cfr$1 = "𝔠";
var Cfr$1 = "ℭ";
var CHcy$1 = "Ч";
var chcy$1 = "ч";
var check$1 = "✓";
var checkmark$1 = "✓";
var Chi$1 = "Χ";
var chi$1 = "χ";
var circ$1 = "ˆ";
var circeq$1 = "≗";
var circlearrowleft$1 = "↺";
var circlearrowright$1 = "↻";
var circledast$1 = "⊛";
var circledcirc$1 = "⊚";
var circleddash$1 = "⊝";
var CircleDot$1 = "⊙";
var circledR$1 = "®";
var circledS$1 = "Ⓢ";
var CircleMinus$1 = "⊖";
var CirclePlus$1 = "⊕";
var CircleTimes$1 = "⊗";
var cir$1 = "○";
var cirE$1 = "⧃";
var cire$1 = "≗";
var cirfnint$1 = "⨐";
var cirmid$1 = "⫯";
var cirscir$1 = "⧂";
var ClockwiseContourIntegral$1 = "∲";
var CloseCurlyDoubleQuote$1 = "”";
var CloseCurlyQuote$1 = "’";
var clubs$1 = "♣";
var clubsuit$1 = "♣";
var colon$1 = ":";
var Colon$1 = "∷";
var Colone$1 = "⩴";
var colone$1 = "≔";
var coloneq$1 = "≔";
var comma$2 = ",";
var commat$1 = "@";
var comp$1 = "∁";
var compfn$1 = "∘";
var complement$1 = "∁";
var complexes$1 = "ℂ";
var cong$1 = "≅";
var congdot$1 = "⩭";
var Congruent$1 = "≡";
var conint$1 = "∮";
var Conint$1 = "∯";
var ContourIntegral$1 = "∮";
var copf$1 = "𝕔";
var Copf$1 = "ℂ";
var coprod$1 = "∐";
var Coproduct$1 = "∐";
var copy$3 = "©";
var COPY$3 = "©";
var copysr$1 = "℗";
var CounterClockwiseContourIntegral$1 = "∳";
var crarr$1 = "↵";
var cross$1 = "✗";
var Cross$1 = "⨯";
var Cscr$1 = "𝒞";
var cscr$1 = "𝒸";
var csub$1 = "⫏";
var csube$1 = "⫑";
var csup$1 = "⫐";
var csupe$1 = "⫒";
var ctdot$1 = "⋯";
var cudarrl$1 = "⤸";
var cudarrr$1 = "⤵";
var cuepr$1 = "⋞";
var cuesc$1 = "⋟";
var cularr$1 = "↶";
var cularrp$1 = "⤽";
var cupbrcap$1 = "⩈";
var cupcap$1 = "⩆";
var CupCap$1 = "≍";
var cup$1 = "∪";
var Cup$1 = "⋓";
var cupcup$1 = "⩊";
var cupdot$1 = "⊍";
var cupor$1 = "⩅";
var cups$1 = "∪︀";
var curarr$1 = "↷";
var curarrm$1 = "⤼";
var curlyeqprec$1 = "⋞";
var curlyeqsucc$1 = "⋟";
var curlyvee$1 = "⋎";
var curlywedge$1 = "⋏";
var curren$3 = "¤";
var curvearrowleft$1 = "↶";
var curvearrowright$1 = "↷";
var cuvee$1 = "⋎";
var cuwed$1 = "⋏";
var cwconint$1 = "∲";
var cwint$1 = "∱";
var cylcty$1 = "⌭";
var dagger$1 = "†";
var Dagger$1 = "‡";
var daleth$1 = "ℸ";
var darr$1 = "↓";
var Darr$1 = "↡";
var dArr$1 = "⇓";
var dash$1 = "‐";
var Dashv$1 = "⫤";
var dashv$1 = "⊣";
var dbkarow$1 = "⤏";
var dblac$1 = "˝";
var Dcaron$1 = "Ď";
var dcaron$1 = "ď";
var Dcy$1 = "Д";
var dcy$1 = "д";
var ddagger$1 = "‡";
var ddarr$1 = "⇊";
var DD$1 = "ⅅ";
var dd$1 = "ⅆ";
var DDotrahd$1 = "⤑";
var ddotseq$1 = "⩷";
var deg$3 = "°";
var Del$1 = "∇";
var Delta$1 = "Δ";
var delta$1 = "δ";
var demptyv$1 = "⦱";
var dfisht$1 = "⥿";
var Dfr$1 = "𝔇";
var dfr$1 = "𝔡";
var dHar$1 = "⥥";
var dharl$1 = "⇃";
var dharr$1 = "⇂";
var DiacriticalAcute$1 = "´";
var DiacriticalDot$1 = "˙";
var DiacriticalDoubleAcute$1 = "˝";
var DiacriticalGrave$1 = "`";
var DiacriticalTilde$1 = "˜";
var diam$1 = "⋄";
var diamond$1 = "⋄";
var Diamond$1 = "⋄";
var diamondsuit$1 = "♦";
var diams$1 = "♦";
var die$1 = "¨";
var DifferentialD$1 = "ⅆ";
var digamma$1 = "ϝ";
var disin$1 = "⋲";
var div$1 = "÷";
var divide$3 = "÷";
var divideontimes$1 = "⋇";
var divonx$1 = "⋇";
var DJcy$1 = "Ђ";
var djcy$1 = "ђ";
var dlcorn$1 = "⌞";
var dlcrop$1 = "⌍";
var dollar$2 = "$";
var Dopf$1 = "𝔻";
var dopf$1 = "𝕕";
var Dot$1 = "¨";
var dot$1 = "˙";
var DotDot$1 = "⃜";
var doteq$1 = "≐";
var doteqdot$1 = "≑";
var DotEqual$1 = "≐";
var dotminus$1 = "∸";
var dotplus$1 = "∔";
var dotsquare$1 = "⊡";
var doublebarwedge$1 = "⌆";
var DoubleContourIntegral$1 = "∯";
var DoubleDot$1 = "¨";
var DoubleDownArrow$1 = "⇓";
var DoubleLeftArrow$1 = "⇐";
var DoubleLeftRightArrow$1 = "⇔";
var DoubleLeftTee$1 = "⫤";
var DoubleLongLeftArrow$1 = "⟸";
var DoubleLongLeftRightArrow$1 = "⟺";
var DoubleLongRightArrow$1 = "⟹";
var DoubleRightArrow$1 = "⇒";
var DoubleRightTee$1 = "⊨";
var DoubleUpArrow$1 = "⇑";
var DoubleUpDownArrow$1 = "⇕";
var DoubleVerticalBar$1 = "∥";
var DownArrowBar$1 = "⤓";
var downarrow$1 = "↓";
var DownArrow$1 = "↓";
var Downarrow$1 = "⇓";
var DownArrowUpArrow$1 = "⇵";
var DownBreve$1 = "̑";
var downdownarrows$1 = "⇊";
var downharpoonleft$1 = "⇃";
var downharpoonright$1 = "⇂";
var DownLeftRightVector$1 = "⥐";
var DownLeftTeeVector$1 = "⥞";
var DownLeftVectorBar$1 = "⥖";
var DownLeftVector$1 = "↽";
var DownRightTeeVector$1 = "⥟";
var DownRightVectorBar$1 = "⥗";
var DownRightVector$1 = "⇁";
var DownTeeArrow$1 = "↧";
var DownTee$1 = "⊤";
var drbkarow$1 = "⤐";
var drcorn$1 = "⌟";
var drcrop$1 = "⌌";
var Dscr$1 = "𝒟";
var dscr$1 = "𝒹";
var DScy$1 = "Ѕ";
var dscy$1 = "ѕ";
var dsol$1 = "⧶";
var Dstrok$1 = "Đ";
var dstrok$1 = "đ";
var dtdot$1 = "⋱";
var dtri$1 = "▿";
var dtrif$1 = "▾";
var duarr$1 = "⇵";
var duhar$1 = "⥯";
var dwangle$1 = "⦦";
var DZcy$1 = "Џ";
var dzcy$1 = "џ";
var dzigrarr$1 = "⟿";
var Eacute$3 = "É";
var eacute$3 = "é";
var easter$1 = "⩮";
var Ecaron$1 = "Ě";
var ecaron$1 = "ě";
var Ecirc$3 = "Ê";
var ecirc$3 = "ê";
var ecir$1 = "≖";
var ecolon$1 = "≕";
var Ecy$1 = "Э";
var ecy$1 = "э";
var eDDot$1 = "⩷";
var Edot$1 = "Ė";
var edot$1 = "ė";
var eDot$1 = "≑";
var ee$1 = "ⅇ";
var efDot$1 = "≒";
var Efr$1 = "𝔈";
var efr$1 = "𝔢";
var eg$1 = "⪚";
var Egrave$3 = "È";
var egrave$3 = "è";
var egs$1 = "⪖";
var egsdot$1 = "⪘";
var el$1 = "⪙";
var Element$1 = "∈";
var elinters$1 = "⏧";
var ell$1 = "ℓ";
var els$1 = "⪕";
var elsdot$1 = "⪗";
var Emacr$1 = "Ē";
var emacr$1 = "ē";
var empty$1 = "∅";
var emptyset$1 = "∅";
var EmptySmallSquare$1 = "◻";
var emptyv$1 = "∅";
var EmptyVerySmallSquare$1 = "▫";
var emsp13$1 = " ";
var emsp14$1 = " ";
var emsp$1 = " ";
var ENG$1 = "Ŋ";
var eng$1 = "ŋ";
var ensp$1 = " ";
var Eogon$1 = "Ę";
var eogon$1 = "ę";
var Eopf$1 = "𝔼";
var eopf$1 = "𝕖";
var epar$1 = "⋕";
var eparsl$1 = "⧣";
var eplus$1 = "⩱";
var epsi$1 = "ε";
var Epsilon$1 = "Ε";
var epsilon$1 = "ε";
var epsiv$1 = "ϵ";
var eqcirc$1 = "≖";
var eqcolon$1 = "≕";
var eqsim$1 = "≂";
var eqslantgtr$1 = "⪖";
var eqslantless$1 = "⪕";
var Equal$1 = "⩵";
var equals$1 = "=";
var EqualTilde$1 = "≂";
var equest$1 = "≟";
var Equilibrium$1 = "⇌";
var equiv$1 = "≡";
var equivDD$1 = "⩸";
var eqvparsl$1 = "⧥";
var erarr$1 = "⥱";
var erDot$1 = "≓";
var escr$1 = "ℯ";
var Escr$1 = "ℰ";
var esdot$1 = "≐";
var Esim$1 = "⩳";
var esim$1 = "≂";
var Eta$1 = "Η";
var eta$1 = "η";
var ETH$3 = "Ð";
var eth$3 = "ð";
var Euml$3 = "Ë";
var euml$3 = "ë";
var euro$1 = "€";
var excl$1 = "!";
var exist$1 = "∃";
var Exists$1 = "∃";
var expectation$1 = "ℰ";
var exponentiale$1 = "ⅇ";
var ExponentialE$1 = "ⅇ";
var fallingdotseq$1 = "≒";
var Fcy$1 = "Ф";
var fcy$1 = "ф";
var female$1 = "♀";
var ffilig$1 = "ﬃ";
var fflig$1 = "ﬀ";
var ffllig$1 = "ﬄ";
var Ffr$1 = "𝔉";
var ffr$1 = "𝔣";
var filig$1 = "ﬁ";
var FilledSmallSquare$1 = "◼";
var FilledVerySmallSquare$1 = "▪";
var fjlig$1 = "fj";
var flat$1 = "♭";
var fllig$1 = "ﬂ";
var fltns$1 = "▱";
var fnof$1 = "ƒ";
var Fopf$1 = "𝔽";
var fopf$1 = "𝕗";
var forall$1 = "∀";
var ForAll$1 = "∀";
var fork$1 = "⋔";
var forkv$1 = "⫙";
var Fouriertrf$1 = "ℱ";
var fpartint$1 = "⨍";
var frac12$3 = "½";
var frac13$1 = "⅓";
var frac14$3 = "¼";
var frac15$1 = "⅕";
var frac16$1 = "⅙";
var frac18$1 = "⅛";
var frac23$1 = "⅔";
var frac25$1 = "⅖";
var frac34$3 = "¾";
var frac35$1 = "⅗";
var frac38$1 = "⅜";
var frac45$1 = "⅘";
var frac56$1 = "⅚";
var frac58$1 = "⅝";
var frac78$1 = "⅞";
var frasl$1 = "⁄";
var frown$1 = "⌢";
var fscr$1 = "𝒻";
var Fscr$1 = "ℱ";
var gacute$1 = "ǵ";
var Gamma$1 = "Γ";
var gamma$1 = "γ";
var Gammad$1 = "Ϝ";
var gammad$1 = "ϝ";
var gap$1 = "⪆";
var Gbreve$1 = "Ğ";
var gbreve$1 = "ğ";
var Gcedil$1 = "Ģ";
var Gcirc$1 = "Ĝ";
var gcirc$1 = "ĝ";
var Gcy$1 = "Г";
var gcy$1 = "г";
var Gdot$1 = "Ġ";
var gdot$1 = "ġ";
var ge$1 = "≥";
var gE$1 = "≧";
var gEl$1 = "⪌";
var gel$1 = "⋛";
var geq$1 = "≥";
var geqq$1 = "≧";
var geqslant$1 = "⩾";
var gescc$1 = "⪩";
var ges$1 = "⩾";
var gesdot$1 = "⪀";
var gesdoto$1 = "⪂";
var gesdotol$1 = "⪄";
var gesl$1 = "⋛︀";
var gesles$1 = "⪔";
var Gfr$1 = "𝔊";
var gfr$1 = "𝔤";
var gg$1 = "≫";
var Gg$1 = "⋙";
var ggg$1 = "⋙";
var gimel$1 = "ℷ";
var GJcy$1 = "Ѓ";
var gjcy$1 = "ѓ";
var gla$1 = "⪥";
var gl$1 = "≷";
var glE$1 = "⪒";
var glj$1 = "⪤";
var gnap$1 = "⪊";
var gnapprox$1 = "⪊";
var gne$1 = "⪈";
var gnE$1 = "≩";
var gneq$1 = "⪈";
var gneqq$1 = "≩";
var gnsim$1 = "⋧";
var Gopf$1 = "𝔾";
var gopf$1 = "𝕘";
var grave$1 = "`";
var GreaterEqual$1 = "≥";
var GreaterEqualLess$1 = "⋛";
var GreaterFullEqual$1 = "≧";
var GreaterGreater$1 = "⪢";
var GreaterLess$1 = "≷";
var GreaterSlantEqual$1 = "⩾";
var GreaterTilde$1 = "≳";
var Gscr$1 = "𝒢";
var gscr$1 = "ℊ";
var gsim$1 = "≳";
var gsime$1 = "⪎";
var gsiml$1 = "⪐";
var gtcc$1 = "⪧";
var gtcir$1 = "⩺";
var gt$6 = ">";
var GT$3 = ">";
var Gt$1 = "≫";
var gtdot$1 = "⋗";
var gtlPar$1 = "⦕";
var gtquest$1 = "⩼";
var gtrapprox$1 = "⪆";
var gtrarr$1 = "⥸";
var gtrdot$1 = "⋗";
var gtreqless$1 = "⋛";
var gtreqqless$1 = "⪌";
var gtrless$1 = "≷";
var gtrsim$1 = "≳";
var gvertneqq$1 = "≩︀";
var gvnE$1 = "≩︀";
var Hacek$1 = "ˇ";
var hairsp$1 = " ";
var half$1 = "½";
var hamilt$1 = "ℋ";
var HARDcy$1 = "Ъ";
var hardcy$1 = "ъ";
var harrcir$1 = "⥈";
var harr$1 = "↔";
var hArr$1 = "⇔";
var harrw$1 = "↭";
var Hat$1 = "^";
var hbar$1 = "ℏ";
var Hcirc$1 = "Ĥ";
var hcirc$1 = "ĥ";
var hearts$1 = "♥";
var heartsuit$1 = "♥";
var hellip$1 = "…";
var hercon$1 = "⊹";
var hfr$1 = "𝔥";
var Hfr$1 = "ℌ";
var HilbertSpace$1 = "ℋ";
var hksearow$1 = "⤥";
var hkswarow$1 = "⤦";
var hoarr$1 = "⇿";
var homtht$1 = "∻";
var hookleftarrow$1 = "↩";
var hookrightarrow$1 = "↪";
var hopf$1 = "𝕙";
var Hopf$1 = "ℍ";
var horbar$1 = "―";
var HorizontalLine$1 = "─";
var hscr$1 = "𝒽";
var Hscr$1 = "ℋ";
var hslash$1 = "ℏ";
var Hstrok$1 = "Ħ";
var hstrok$1 = "ħ";
var HumpDownHump$1 = "≎";
var HumpEqual$1 = "≏";
var hybull$1 = "⁃";
var hyphen$1 = "‐";
var Iacute$3 = "Í";
var iacute$3 = "í";
var ic$1 = "⁣";
var Icirc$3 = "Î";
var icirc$3 = "î";
var Icy$1 = "И";
var icy$1 = "и";
var Idot$1 = "İ";
var IEcy$1 = "Е";
var iecy$1 = "е";
var iexcl$3 = "¡";
var iff$1 = "⇔";
var ifr$1 = "𝔦";
var Ifr$1 = "ℑ";
var Igrave$3 = "Ì";
var igrave$3 = "ì";
var ii$1 = "ⅈ";
var iiiint$1 = "⨌";
var iiint$1 = "∭";
var iinfin$1 = "⧜";
var iiota$1 = "℩";
var IJlig$1 = "Ĳ";
var ijlig$1 = "ĳ";
var Imacr$1 = "Ī";
var imacr$1 = "ī";
var image$1 = "ℑ";
var ImaginaryI$1 = "ⅈ";
var imagline$1 = "ℐ";
var imagpart$1 = "ℑ";
var imath$1 = "ı";
var Im$1 = "ℑ";
var imof$1 = "⊷";
var imped$1 = "Ƶ";
var Implies$1 = "⇒";
var incare$1 = "℅";
var infin$1 = "∞";
var infintie$1 = "⧝";
var inodot$1 = "ı";
var intcal$1 = "⊺";
var int$1 = "∫";
var Int$1 = "∬";
var integers$1 = "ℤ";
var Integral$1 = "∫";
var intercal$1 = "⊺";
var Intersection$1 = "⋂";
var intlarhk$1 = "⨗";
var intprod$1 = "⨼";
var InvisibleComma$1 = "⁣";
var InvisibleTimes$1 = "⁢";
var IOcy$1 = "Ё";
var iocy$1 = "ё";
var Iogon$1 = "Į";
var iogon$1 = "į";
var Iopf$1 = "𝕀";
var iopf$1 = "𝕚";
var Iota$1 = "Ι";
var iota$1 = "ι";
var iprod$1 = "⨼";
var iquest$3 = "¿";
var iscr$1 = "𝒾";
var Iscr$1 = "ℐ";
var isin$1 = "∈";
var isindot$1 = "⋵";
var isinE$1 = "⋹";
var isins$1 = "⋴";
var isinsv$1 = "⋳";
var isinv$1 = "∈";
var it$1 = "⁢";
var Itilde$1 = "Ĩ";
var itilde$1 = "ĩ";
var Iukcy$1 = "І";
var iukcy$1 = "і";
var Iuml$3 = "Ï";
var iuml$3 = "ï";
var Jcirc$1 = "Ĵ";
var jcirc$1 = "ĵ";
var Jcy$1 = "Й";
var jcy$1 = "й";
var Jfr$1 = "𝔍";
var jfr$1 = "𝔧";
var jmath$1 = "ȷ";
var Jopf$1 = "𝕁";
var jopf$1 = "𝕛";
var Jscr$1 = "𝒥";
var jscr$1 = "𝒿";
var Jsercy$1 = "Ј";
var jsercy$1 = "ј";
var Jukcy$1 = "Є";
var jukcy$1 = "є";
var Kappa$1 = "Κ";
var kappa$1 = "κ";
var kappav$1 = "ϰ";
var Kcedil$1 = "Ķ";
var kcedil$1 = "ķ";
var Kcy$1 = "К";
var kcy$1 = "к";
var Kfr$1 = "𝔎";
var kfr$1 = "𝔨";
var kgreen$1 = "ĸ";
var KHcy$1 = "Х";
var khcy$1 = "х";
var KJcy$1 = "Ќ";
var kjcy$1 = "ќ";
var Kopf$1 = "𝕂";
var kopf$1 = "𝕜";
var Kscr$1 = "𝒦";
var kscr$1 = "𝓀";
var lAarr$1 = "⇚";
var Lacute$1 = "Ĺ";
var lacute$1 = "ĺ";
var laemptyv$1 = "⦴";
var lagran$1 = "ℒ";
var Lambda$1 = "Λ";
var lambda$1 = "λ";
var lang$1 = "⟨";
var Lang$1 = "⟪";
var langd$1 = "⦑";
var langle$1 = "⟨";
var lap$1 = "⪅";
var Laplacetrf$1 = "ℒ";
var laquo$3 = "«";
var larrb$1 = "⇤";
var larrbfs$1 = "⤟";
var larr$1 = "←";
var Larr$1 = "↞";
var lArr$1 = "⇐";
var larrfs$1 = "⤝";
var larrhk$1 = "↩";
var larrlp$1 = "↫";
var larrpl$1 = "⤹";
var larrsim$1 = "⥳";
var larrtl$1 = "↢";
var latail$1 = "⤙";
var lAtail$1 = "⤛";
var lat$1 = "⪫";
var late$1 = "⪭";
var lates$1 = "⪭︀";
var lbarr$1 = "⤌";
var lBarr$1 = "⤎";
var lbbrk$1 = "❲";
var lbrace$1 = "{";
var lbrack$1 = "[";
var lbrke$1 = "⦋";
var lbrksld$1 = "⦏";
var lbrkslu$1 = "⦍";
var Lcaron$1 = "Ľ";
var lcaron$1 = "ľ";
var Lcedil$1 = "Ļ";
var lcedil$1 = "ļ";
var lceil$1 = "⌈";
var lcub$1 = "{";
var Lcy$1 = "Л";
var lcy$1 = "л";
var ldca$1 = "⤶";
var ldquo$1 = "“";
var ldquor$1 = "„";
var ldrdhar$1 = "⥧";
var ldrushar$1 = "⥋";
var ldsh$1 = "↲";
var le$1 = "≤";
var lE$1 = "≦";
var LeftAngleBracket$1 = "⟨";
var LeftArrowBar$1 = "⇤";
var leftarrow$1 = "←";
var LeftArrow$1 = "←";
var Leftarrow$1 = "⇐";
var LeftArrowRightArrow$1 = "⇆";
var leftarrowtail$1 = "↢";
var LeftCeiling$1 = "⌈";
var LeftDoubleBracket$1 = "⟦";
var LeftDownTeeVector$1 = "⥡";
var LeftDownVectorBar$1 = "⥙";
var LeftDownVector$1 = "⇃";
var LeftFloor$1 = "⌊";
var leftharpoondown$1 = "↽";
var leftharpoonup$1 = "↼";
var leftleftarrows$1 = "⇇";
var leftrightarrow$1 = "↔";
var LeftRightArrow$1 = "↔";
var Leftrightarrow$1 = "⇔";
var leftrightarrows$1 = "⇆";
var leftrightharpoons$1 = "⇋";
var leftrightsquigarrow$1 = "↭";
var LeftRightVector$1 = "⥎";
var LeftTeeArrow$1 = "↤";
var LeftTee$1 = "⊣";
var LeftTeeVector$1 = "⥚";
var leftthreetimes$1 = "⋋";
var LeftTriangleBar$1 = "⧏";
var LeftTriangle$1 = "⊲";
var LeftTriangleEqual$1 = "⊴";
var LeftUpDownVector$1 = "⥑";
var LeftUpTeeVector$1 = "⥠";
var LeftUpVectorBar$1 = "⥘";
var LeftUpVector$1 = "↿";
var LeftVectorBar$1 = "⥒";
var LeftVector$1 = "↼";
var lEg$1 = "⪋";
var leg$1 = "⋚";
var leq$1 = "≤";
var leqq$1 = "≦";
var leqslant$1 = "⩽";
var lescc$1 = "⪨";
var les$1 = "⩽";
var lesdot$1 = "⩿";
var lesdoto$1 = "⪁";
var lesdotor$1 = "⪃";
var lesg$1 = "⋚︀";
var lesges$1 = "⪓";
var lessapprox$1 = "⪅";
var lessdot$1 = "⋖";
var lesseqgtr$1 = "⋚";
var lesseqqgtr$1 = "⪋";
var LessEqualGreater$1 = "⋚";
var LessFullEqual$1 = "≦";
var LessGreater$1 = "≶";
var lessgtr$1 = "≶";
var LessLess$1 = "⪡";
var lesssim$1 = "≲";
var LessSlantEqual$1 = "⩽";
var LessTilde$1 = "≲";
var lfisht$1 = "⥼";
var lfloor$1 = "⌊";
var Lfr$1 = "𝔏";
var lfr$1 = "𝔩";
var lg$1 = "≶";
var lgE$1 = "⪑";
var lHar$1 = "⥢";
var lhard$1 = "↽";
var lharu$1 = "↼";
var lharul$1 = "⥪";
var lhblk$1 = "▄";
var LJcy$1 = "Љ";
var ljcy$1 = "љ";
var llarr$1 = "⇇";
var ll$1 = "≪";
var Ll$1 = "⋘";
var llcorner$1 = "⌞";
var Lleftarrow$1 = "⇚";
var llhard$1 = "⥫";
var lltri$1 = "◺";
var Lmidot$1 = "Ŀ";
var lmidot$1 = "ŀ";
var lmoustache$1 = "⎰";
var lmoust$1 = "⎰";
var lnap$1 = "⪉";
var lnapprox$1 = "⪉";
var lne$1 = "⪇";
var lnE$1 = "≨";
var lneq$1 = "⪇";
var lneqq$1 = "≨";
var lnsim$1 = "⋦";
var loang$1 = "⟬";
var loarr$1 = "⇽";
var lobrk$1 = "⟦";
var longleftarrow$1 = "⟵";
var LongLeftArrow$1 = "⟵";
var Longleftarrow$1 = "⟸";
var longleftrightarrow$1 = "⟷";
var LongLeftRightArrow$1 = "⟷";
var Longleftrightarrow$1 = "⟺";
var longmapsto$1 = "⟼";
var longrightarrow$1 = "⟶";
var LongRightArrow$1 = "⟶";
var Longrightarrow$1 = "⟹";
var looparrowleft$1 = "↫";
var looparrowright$1 = "↬";
var lopar$1 = "⦅";
var Lopf$1 = "𝕃";
var lopf$1 = "𝕝";
var loplus$1 = "⨭";
var lotimes$1 = "⨴";
var lowast$1 = "∗";
var lowbar$1 = "_";
var LowerLeftArrow$1 = "↙";
var LowerRightArrow$1 = "↘";
var loz$1 = "◊";
var lozenge$1 = "◊";
var lozf$1 = "⧫";
var lpar$1 = "(";
var lparlt$1 = "⦓";
var lrarr$1 = "⇆";
var lrcorner$1 = "⌟";
var lrhar$1 = "⇋";
var lrhard$1 = "⥭";
var lrm$1 = "‎";
var lrtri$1 = "⊿";
var lsaquo$1 = "‹";
var lscr$1 = "𝓁";
var Lscr$1 = "ℒ";
var lsh$1 = "↰";
var Lsh$1 = "↰";
var lsim$1 = "≲";
var lsime$1 = "⪍";
var lsimg$1 = "⪏";
var lsqb$1 = "[";
var lsquo$1 = "‘";
var lsquor$1 = "‚";
var Lstrok$1 = "Ł";
var lstrok$1 = "ł";
var ltcc$1 = "⪦";
var ltcir$1 = "⩹";
var lt$5 = "<";
var LT$3 = "<";
var Lt$1 = "≪";
var ltdot$1 = "⋖";
var lthree$1 = "⋋";
var ltimes$1 = "⋉";
var ltlarr$1 = "⥶";
var ltquest$1 = "⩻";
var ltri$1 = "◃";
var ltrie$1 = "⊴";
var ltrif$1 = "◂";
var ltrPar$1 = "⦖";
var lurdshar$1 = "⥊";
var luruhar$1 = "⥦";
var lvertneqq$1 = "≨︀";
var lvnE$1 = "≨︀";
var macr$3 = "¯";
var male$1 = "♂";
var malt$1 = "✠";
var maltese$1 = "✠";
var map$1 = "↦";
var mapsto$1 = "↦";
var mapstodown$1 = "↧";
var mapstoleft$1 = "↤";
var mapstoup$1 = "↥";
var marker$1 = "▮";
var mcomma$1 = "⨩";
var Mcy$1 = "М";
var mcy$1 = "м";
var mdash$1 = "—";
var mDDot$1 = "∺";
var measuredangle$1 = "∡";
var MediumSpace$1 = " ";
var Mellintrf$1 = "ℳ";
var Mfr$1 = "𝔐";
var mfr$1 = "𝔪";
var mho$1 = "℧";
var micro$3 = "µ";
var midast$1 = "*";
var midcir$1 = "⫰";
var mid$1 = "∣";
var middot$3 = "·";
var minusb$1 = "⊟";
var minus$1 = "−";
var minusd$1 = "∸";
var minusdu$1 = "⨪";
var MinusPlus$1 = "∓";
var mlcp$1 = "⫛";
var mldr$1 = "…";
var mnplus$1 = "∓";
var models$1 = "⊧";
var Mopf$1 = "𝕄";
var mopf$1 = "𝕞";
var mp$1 = "∓";
var mscr$1 = "𝓂";
var Mscr$1 = "ℳ";
var mstpos$1 = "∾";
var Mu$1 = "Μ";
var mu$1 = "μ";
var multimap$1 = "⊸";
var mumap$1 = "⊸";
var nabla$1 = "∇";
var Nacute$1 = "Ń";
var nacute$1 = "ń";
var nang$1 = "∠⃒";
var nap$1 = "≉";
var napE$1 = "⩰̸";
var napid$1 = "≋̸";
var napos$1 = "ŉ";
var napprox$1 = "≉";
var natural$1 = "♮";
var naturals$1 = "ℕ";
var natur$1 = "♮";
var nbsp$3 = " ";
var nbump$1 = "≎̸";
var nbumpe$1 = "≏̸";
var ncap$1 = "⩃";
var Ncaron$1 = "Ň";
var ncaron$1 = "ň";
var Ncedil$1 = "Ņ";
var ncedil$1 = "ņ";
var ncong$1 = "≇";
var ncongdot$1 = "⩭̸";
var ncup$1 = "⩂";
var Ncy$1 = "Н";
var ncy$1 = "н";
var ndash$1 = "–";
var nearhk$1 = "⤤";
var nearr$1 = "↗";
var neArr$1 = "⇗";
var nearrow$1 = "↗";
var ne$1 = "≠";
var nedot$1 = "≐̸";
var NegativeMediumSpace$1 = "​";
var NegativeThickSpace$1 = "​";
var NegativeThinSpace$1 = "​";
var NegativeVeryThinSpace$1 = "​";
var nequiv$1 = "≢";
var nesear$1 = "⤨";
var nesim$1 = "≂̸";
var NestedGreaterGreater$1 = "≫";
var NestedLessLess$1 = "≪";
var NewLine$1 = "\n";
var nexist$1 = "∄";
var nexists$1 = "∄";
var Nfr$1 = "𝔑";
var nfr$1 = "𝔫";
var ngE$1 = "≧̸";
var nge$1 = "≱";
var ngeq$1 = "≱";
var ngeqq$1 = "≧̸";
var ngeqslant$1 = "⩾̸";
var nges$1 = "⩾̸";
var nGg$1 = "⋙̸";
var ngsim$1 = "≵";
var nGt$1 = "≫⃒";
var ngt$1 = "≯";
var ngtr$1 = "≯";
var nGtv$1 = "≫̸";
var nharr$1 = "↮";
var nhArr$1 = "⇎";
var nhpar$1 = "⫲";
var ni$1 = "∋";
var nis$1 = "⋼";
var nisd$1 = "⋺";
var niv$1 = "∋";
var NJcy$1 = "Њ";
var njcy$1 = "њ";
var nlarr$1 = "↚";
var nlArr$1 = "⇍";
var nldr$1 = "‥";
var nlE$1 = "≦̸";
var nle$1 = "≰";
var nleftarrow$1 = "↚";
var nLeftarrow$1 = "⇍";
var nleftrightarrow$1 = "↮";
var nLeftrightarrow$1 = "⇎";
var nleq$1 = "≰";
var nleqq$1 = "≦̸";
var nleqslant$1 = "⩽̸";
var nles$1 = "⩽̸";
var nless$1 = "≮";
var nLl$1 = "⋘̸";
var nlsim$1 = "≴";
var nLt$1 = "≪⃒";
var nlt$1 = "≮";
var nltri$1 = "⋪";
var nltrie$1 = "⋬";
var nLtv$1 = "≪̸";
var nmid$1 = "∤";
var NoBreak$1 = "⁠";
var NonBreakingSpace$1 = " ";
var nopf$1 = "𝕟";
var Nopf$1 = "ℕ";
var Not$1 = "⫬";
var not$3 = "¬";
var NotCongruent$1 = "≢";
var NotCupCap$1 = "≭";
var NotDoubleVerticalBar$1 = "∦";
var NotElement$1 = "∉";
var NotEqual$1 = "≠";
var NotEqualTilde$1 = "≂̸";
var NotExists$1 = "∄";
var NotGreater$1 = "≯";
var NotGreaterEqual$1 = "≱";
var NotGreaterFullEqual$1 = "≧̸";
var NotGreaterGreater$1 = "≫̸";
var NotGreaterLess$1 = "≹";
var NotGreaterSlantEqual$1 = "⩾̸";
var NotGreaterTilde$1 = "≵";
var NotHumpDownHump$1 = "≎̸";
var NotHumpEqual$1 = "≏̸";
var notin$1 = "∉";
var notindot$1 = "⋵̸";
var notinE$1 = "⋹̸";
var notinva$1 = "∉";
var notinvb$1 = "⋷";
var notinvc$1 = "⋶";
var NotLeftTriangleBar$1 = "⧏̸";
var NotLeftTriangle$1 = "⋪";
var NotLeftTriangleEqual$1 = "⋬";
var NotLess$1 = "≮";
var NotLessEqual$1 = "≰";
var NotLessGreater$1 = "≸";
var NotLessLess$1 = "≪̸";
var NotLessSlantEqual$1 = "⩽̸";
var NotLessTilde$1 = "≴";
var NotNestedGreaterGreater$1 = "⪢̸";
var NotNestedLessLess$1 = "⪡̸";
var notni$1 = "∌";
var notniva$1 = "∌";
var notnivb$1 = "⋾";
var notnivc$1 = "⋽";
var NotPrecedes$1 = "⊀";
var NotPrecedesEqual$1 = "⪯̸";
var NotPrecedesSlantEqual$1 = "⋠";
var NotReverseElement$1 = "∌";
var NotRightTriangleBar$1 = "⧐̸";
var NotRightTriangle$1 = "⋫";
var NotRightTriangleEqual$1 = "⋭";
var NotSquareSubset$1 = "⊏̸";
var NotSquareSubsetEqual$1 = "⋢";
var NotSquareSuperset$1 = "⊐̸";
var NotSquareSupersetEqual$1 = "⋣";
var NotSubset$1 = "⊂⃒";
var NotSubsetEqual$1 = "⊈";
var NotSucceeds$1 = "⊁";
var NotSucceedsEqual$1 = "⪰̸";
var NotSucceedsSlantEqual$1 = "⋡";
var NotSucceedsTilde$1 = "≿̸";
var NotSuperset$1 = "⊃⃒";
var NotSupersetEqual$1 = "⊉";
var NotTilde$1 = "≁";
var NotTildeEqual$1 = "≄";
var NotTildeFullEqual$1 = "≇";
var NotTildeTilde$1 = "≉";
var NotVerticalBar$1 = "∤";
var nparallel$1 = "∦";
var npar$1 = "∦";
var nparsl$1 = "⫽⃥";
var npart$1 = "∂̸";
var npolint$1 = "⨔";
var npr$1 = "⊀";
var nprcue$1 = "⋠";
var nprec$1 = "⊀";
var npreceq$1 = "⪯̸";
var npre$1 = "⪯̸";
var nrarrc$1 = "⤳̸";
var nrarr$1 = "↛";
var nrArr$1 = "⇏";
var nrarrw$1 = "↝̸";
var nrightarrow$1 = "↛";
var nRightarrow$1 = "⇏";
var nrtri$1 = "⋫";
var nrtrie$1 = "⋭";
var nsc$1 = "⊁";
var nsccue$1 = "⋡";
var nsce$1 = "⪰̸";
var Nscr$1 = "𝒩";
var nscr$1 = "𝓃";
var nshortmid$1 = "∤";
var nshortparallel$1 = "∦";
var nsim$1 = "≁";
var nsime$1 = "≄";
var nsimeq$1 = "≄";
var nsmid$1 = "∤";
var nspar$1 = "∦";
var nsqsube$1 = "⋢";
var nsqsupe$1 = "⋣";
var nsub$1 = "⊄";
var nsubE$1 = "⫅̸";
var nsube$1 = "⊈";
var nsubset$1 = "⊂⃒";
var nsubseteq$1 = "⊈";
var nsubseteqq$1 = "⫅̸";
var nsucc$1 = "⊁";
var nsucceq$1 = "⪰̸";
var nsup$1 = "⊅";
var nsupE$1 = "⫆̸";
var nsupe$1 = "⊉";
var nsupset$1 = "⊃⃒";
var nsupseteq$1 = "⊉";
var nsupseteqq$1 = "⫆̸";
var ntgl$1 = "≹";
var Ntilde$3 = "Ñ";
var ntilde$3 = "ñ";
var ntlg$1 = "≸";
var ntriangleleft$1 = "⋪";
var ntrianglelefteq$1 = "⋬";
var ntriangleright$1 = "⋫";
var ntrianglerighteq$1 = "⋭";
var Nu$1 = "Ν";
var nu$1 = "ν";
var num$1 = "#";
var numero$1 = "№";
var numsp$1 = " ";
var nvap$1 = "≍⃒";
var nvdash$1 = "⊬";
var nvDash$1 = "⊭";
var nVdash$1 = "⊮";
var nVDash$1 = "⊯";
var nvge$1 = "≥⃒";
var nvgt$1 = ">⃒";
var nvHarr$1 = "⤄";
var nvinfin$1 = "⧞";
var nvlArr$1 = "⤂";
var nvle$1 = "≤⃒";
var nvlt$1 = "<⃒";
var nvltrie$1 = "⊴⃒";
var nvrArr$1 = "⤃";
var nvrtrie$1 = "⊵⃒";
var nvsim$1 = "∼⃒";
var nwarhk$1 = "⤣";
var nwarr$1 = "↖";
var nwArr$1 = "⇖";
var nwarrow$1 = "↖";
var nwnear$1 = "⤧";
var Oacute$3 = "Ó";
var oacute$3 = "ó";
var oast$1 = "⊛";
var Ocirc$3 = "Ô";
var ocirc$3 = "ô";
var ocir$1 = "⊚";
var Ocy$1 = "О";
var ocy$1 = "о";
var odash$1 = "⊝";
var Odblac$1 = "Ő";
var odblac$1 = "ő";
var odiv$1 = "⨸";
var odot$1 = "⊙";
var odsold$1 = "⦼";
var OElig$1 = "Œ";
var oelig$1 = "œ";
var ofcir$1 = "⦿";
var Ofr$1 = "𝔒";
var ofr$1 = "𝔬";
var ogon$1 = "˛";
var Ograve$3 = "Ò";
var ograve$3 = "ò";
var ogt$1 = "⧁";
var ohbar$1 = "⦵";
var ohm$1 = "Ω";
var oint$1 = "∮";
var olarr$1 = "↺";
var olcir$1 = "⦾";
var olcross$1 = "⦻";
var oline$1 = "‾";
var olt$1 = "⧀";
var Omacr$1 = "Ō";
var omacr$1 = "ō";
var Omega$1 = "Ω";
var omega$1 = "ω";
var Omicron$1 = "Ο";
var omicron$1 = "ο";
var omid$1 = "⦶";
var ominus$1 = "⊖";
var Oopf$1 = "𝕆";
var oopf$1 = "𝕠";
var opar$1 = "⦷";
var OpenCurlyDoubleQuote$1 = "“";
var OpenCurlyQuote$1 = "‘";
var operp$1 = "⦹";
var oplus$1 = "⊕";
var orarr$1 = "↻";
var Or$1 = "⩔";
var or$1 = "∨";
var ord$1 = "⩝";
var order$1 = "ℴ";
var orderof$1 = "ℴ";
var ordf$3 = "ª";
var ordm$3 = "º";
var origof$1 = "⊶";
var oror$1 = "⩖";
var orslope$1 = "⩗";
var orv$1 = "⩛";
var oS$1 = "Ⓢ";
var Oscr$1 = "𝒪";
var oscr$1 = "ℴ";
var Oslash$3 = "Ø";
var oslash$3 = "ø";
var osol$1 = "⊘";
var Otilde$3 = "Õ";
var otilde$3 = "õ";
var otimesas$1 = "⨶";
var Otimes$1 = "⨷";
var otimes$1 = "⊗";
var Ouml$3 = "Ö";
var ouml$3 = "ö";
var ovbar$1 = "⌽";
var OverBar$1 = "‾";
var OverBrace$1 = "⏞";
var OverBracket$1 = "⎴";
var OverParenthesis$1 = "⏜";
var para$3 = "¶";
var parallel$1 = "∥";
var par$1 = "∥";
var parsim$1 = "⫳";
var parsl$1 = "⫽";
var part$1 = "∂";
var PartialD$1 = "∂";
var Pcy$1 = "П";
var pcy$1 = "п";
var percnt$1 = "%";
var period$1 = ".";
var permil$1 = "‰";
var perp$1 = "⊥";
var pertenk$1 = "‱";
var Pfr$1 = "𝔓";
var pfr$1 = "𝔭";
var Phi$1 = "Φ";
var phi$1 = "φ";
var phiv$1 = "ϕ";
var phmmat$1 = "ℳ";
var phone$1 = "☎";
var Pi$1 = "Π";
var pi$1 = "π";
var pitchfork$1 = "⋔";
var piv$1 = "ϖ";
var planck$1 = "ℏ";
var planckh$1 = "ℎ";
var plankv$1 = "ℏ";
var plusacir$1 = "⨣";
var plusb$1 = "⊞";
var pluscir$1 = "⨢";
var plus$2 = "+";
var plusdo$1 = "∔";
var plusdu$1 = "⨥";
var pluse$1 = "⩲";
var PlusMinus$1 = "±";
var plusmn$3 = "±";
var plussim$1 = "⨦";
var plustwo$1 = "⨧";
var pm$1 = "±";
var Poincareplane$1 = "ℌ";
var pointint$1 = "⨕";
var popf$1 = "𝕡";
var Popf$1 = "ℙ";
var pound$3 = "£";
var prap$1 = "⪷";
var Pr$1 = "⪻";
var pr$1 = "≺";
var prcue$1 = "≼";
var precapprox$1 = "⪷";
var prec$1 = "≺";
var preccurlyeq$1 = "≼";
var Precedes$1 = "≺";
var PrecedesEqual$1 = "⪯";
var PrecedesSlantEqual$1 = "≼";
var PrecedesTilde$1 = "≾";
var preceq$1 = "⪯";
var precnapprox$1 = "⪹";
var precneqq$1 = "⪵";
var precnsim$1 = "⋨";
var pre$1 = "⪯";
var prE$1 = "⪳";
var precsim$1 = "≾";
var prime$1 = "′";
var Prime$1 = "″";
var primes$1 = "ℙ";
var prnap$1 = "⪹";
var prnE$1 = "⪵";
var prnsim$1 = "⋨";
var prod$1 = "∏";
var Product$1 = "∏";
var profalar$1 = "⌮";
var profline$1 = "⌒";
var profsurf$1 = "⌓";
var prop$1 = "∝";
var Proportional$1 = "∝";
var Proportion$1 = "∷";
var propto$1 = "∝";
var prsim$1 = "≾";
var prurel$1 = "⊰";
var Pscr$1 = "𝒫";
var pscr$1 = "𝓅";
var Psi$1 = "Ψ";
var psi$1 = "ψ";
var puncsp$1 = " ";
var Qfr$1 = "𝔔";
var qfr$1 = "𝔮";
var qint$1 = "⨌";
var qopf$1 = "𝕢";
var Qopf$1 = "ℚ";
var qprime$1 = "⁗";
var Qscr$1 = "𝒬";
var qscr$1 = "𝓆";
var quaternions$1 = "ℍ";
var quatint$1 = "⨖";
var quest$1 = "?";
var questeq$1 = "≟";
var quot$5 = "\"";
var QUOT$3 = "\"";
var rAarr$1 = "⇛";
var race$1 = "∽̱";
var Racute$1 = "Ŕ";
var racute$1 = "ŕ";
var radic$1 = "√";
var raemptyv$1 = "⦳";
var rang$1 = "⟩";
var Rang$1 = "⟫";
var rangd$1 = "⦒";
var range$1 = "⦥";
var rangle$1 = "⟩";
var raquo$3 = "»";
var rarrap$1 = "⥵";
var rarrb$1 = "⇥";
var rarrbfs$1 = "⤠";
var rarrc$1 = "⤳";
var rarr$1 = "→";
var Rarr$1 = "↠";
var rArr$1 = "⇒";
var rarrfs$1 = "⤞";
var rarrhk$1 = "↪";
var rarrlp$1 = "↬";
var rarrpl$1 = "⥅";
var rarrsim$1 = "⥴";
var Rarrtl$1 = "⤖";
var rarrtl$1 = "↣";
var rarrw$1 = "↝";
var ratail$1 = "⤚";
var rAtail$1 = "⤜";
var ratio$1 = "∶";
var rationals$1 = "ℚ";
var rbarr$1 = "⤍";
var rBarr$1 = "⤏";
var RBarr$1 = "⤐";
var rbbrk$1 = "❳";
var rbrace$1 = "}";
var rbrack$1 = "]";
var rbrke$1 = "⦌";
var rbrksld$1 = "⦎";
var rbrkslu$1 = "⦐";
var Rcaron$1 = "Ř";
var rcaron$1 = "ř";
var Rcedil$1 = "Ŗ";
var rcedil$1 = "ŗ";
var rceil$1 = "⌉";
var rcub$1 = "}";
var Rcy$1 = "Р";
var rcy$1 = "р";
var rdca$1 = "⤷";
var rdldhar$1 = "⥩";
var rdquo$1 = "”";
var rdquor$1 = "”";
var rdsh$1 = "↳";
var real$1 = "ℜ";
var realine$1 = "ℛ";
var realpart$1 = "ℜ";
var reals$1 = "ℝ";
var Re$1 = "ℜ";
var rect$1 = "▭";
var reg$3 = "®";
var REG$3 = "®";
var ReverseElement$1 = "∋";
var ReverseEquilibrium$1 = "⇋";
var ReverseUpEquilibrium$1 = "⥯";
var rfisht$1 = "⥽";
var rfloor$1 = "⌋";
var rfr$1 = "𝔯";
var Rfr$1 = "ℜ";
var rHar$1 = "⥤";
var rhard$1 = "⇁";
var rharu$1 = "⇀";
var rharul$1 = "⥬";
var Rho$1 = "Ρ";
var rho$1 = "ρ";
var rhov$1 = "ϱ";
var RightAngleBracket$1 = "⟩";
var RightArrowBar$1 = "⇥";
var rightarrow$1 = "→";
var RightArrow$1 = "→";
var Rightarrow$1 = "⇒";
var RightArrowLeftArrow$1 = "⇄";
var rightarrowtail$1 = "↣";
var RightCeiling$1 = "⌉";
var RightDoubleBracket$1 = "⟧";
var RightDownTeeVector$1 = "⥝";
var RightDownVectorBar$1 = "⥕";
var RightDownVector$1 = "⇂";
var RightFloor$1 = "⌋";
var rightharpoondown$1 = "⇁";
var rightharpoonup$1 = "⇀";
var rightleftarrows$1 = "⇄";
var rightleftharpoons$1 = "⇌";
var rightrightarrows$1 = "⇉";
var rightsquigarrow$1 = "↝";
var RightTeeArrow$1 = "↦";
var RightTee$1 = "⊢";
var RightTeeVector$1 = "⥛";
var rightthreetimes$1 = "⋌";
var RightTriangleBar$1 = "⧐";
var RightTriangle$1 = "⊳";
var RightTriangleEqual$1 = "⊵";
var RightUpDownVector$1 = "⥏";
var RightUpTeeVector$1 = "⥜";
var RightUpVectorBar$1 = "⥔";
var RightUpVector$1 = "↾";
var RightVectorBar$1 = "⥓";
var RightVector$1 = "⇀";
var ring$1 = "˚";
var risingdotseq$1 = "≓";
var rlarr$1 = "⇄";
var rlhar$1 = "⇌";
var rlm$1 = "‏";
var rmoustache$1 = "⎱";
var rmoust$1 = "⎱";
var rnmid$1 = "⫮";
var roang$1 = "⟭";
var roarr$1 = "⇾";
var robrk$1 = "⟧";
var ropar$1 = "⦆";
var ropf$1 = "𝕣";
var Ropf$1 = "ℝ";
var roplus$1 = "⨮";
var rotimes$1 = "⨵";
var RoundImplies$1 = "⥰";
var rpar$1 = ")";
var rpargt$1 = "⦔";
var rppolint$1 = "⨒";
var rrarr$1 = "⇉";
var Rrightarrow$1 = "⇛";
var rsaquo$1 = "›";
var rscr$1 = "𝓇";
var Rscr$1 = "ℛ";
var rsh$1 = "↱";
var Rsh$1 = "↱";
var rsqb$1 = "]";
var rsquo$1 = "’";
var rsquor$1 = "’";
var rthree$1 = "⋌";
var rtimes$1 = "⋊";
var rtri$1 = "▹";
var rtrie$1 = "⊵";
var rtrif$1 = "▸";
var rtriltri$1 = "⧎";
var RuleDelayed$1 = "⧴";
var ruluhar$1 = "⥨";
var rx$1 = "℞";
var Sacute$1 = "Ś";
var sacute$1 = "ś";
var sbquo$1 = "‚";
var scap$1 = "⪸";
var Scaron$1 = "Š";
var scaron$1 = "š";
var Sc$1 = "⪼";
var sc$1 = "≻";
var sccue$1 = "≽";
var sce$1 = "⪰";
var scE$1 = "⪴";
var Scedil$1 = "Ş";
var scedil$1 = "ş";
var Scirc$1 = "Ŝ";
var scirc$1 = "ŝ";
var scnap$1 = "⪺";
var scnE$1 = "⪶";
var scnsim$1 = "⋩";
var scpolint$1 = "⨓";
var scsim$1 = "≿";
var Scy$1 = "С";
var scy$1 = "с";
var sdotb$1 = "⊡";
var sdot$1 = "⋅";
var sdote$1 = "⩦";
var searhk$1 = "⤥";
var searr$1 = "↘";
var seArr$1 = "⇘";
var searrow$1 = "↘";
var sect$3 = "§";
var semi$1 = ";";
var seswar$1 = "⤩";
var setminus$1 = "∖";
var setmn$1 = "∖";
var sext$1 = "✶";
var Sfr$1 = "𝔖";
var sfr$1 = "𝔰";
var sfrown$1 = "⌢";
var sharp$1 = "♯";
var SHCHcy$1 = "Щ";
var shchcy$1 = "щ";
var SHcy$1 = "Ш";
var shcy$1 = "ш";
var ShortDownArrow$1 = "↓";
var ShortLeftArrow$1 = "←";
var shortmid$1 = "∣";
var shortparallel$1 = "∥";
var ShortRightArrow$1 = "→";
var ShortUpArrow$1 = "↑";
var shy$3 = "­";
var Sigma$1 = "Σ";
var sigma$1 = "σ";
var sigmaf$1 = "ς";
var sigmav$1 = "ς";
var sim$1 = "∼";
var simdot$1 = "⩪";
var sime$1 = "≃";
var simeq$1 = "≃";
var simg$1 = "⪞";
var simgE$1 = "⪠";
var siml$1 = "⪝";
var simlE$1 = "⪟";
var simne$1 = "≆";
var simplus$1 = "⨤";
var simrarr$1 = "⥲";
var slarr$1 = "←";
var SmallCircle$1 = "∘";
var smallsetminus$1 = "∖";
var smashp$1 = "⨳";
var smeparsl$1 = "⧤";
var smid$1 = "∣";
var smile$1 = "⌣";
var smt$1 = "⪪";
var smte$1 = "⪬";
var smtes$1 = "⪬︀";
var SOFTcy$1 = "Ь";
var softcy$1 = "ь";
var solbar$1 = "⌿";
var solb$1 = "⧄";
var sol$1 = "/";
var Sopf$1 = "𝕊";
var sopf$1 = "𝕤";
var spades$1 = "♠";
var spadesuit$1 = "♠";
var spar$1 = "∥";
var sqcap$1 = "⊓";
var sqcaps$1 = "⊓︀";
var sqcup$1 = "⊔";
var sqcups$1 = "⊔︀";
var Sqrt$1 = "√";
var sqsub$1 = "⊏";
var sqsube$1 = "⊑";
var sqsubset$1 = "⊏";
var sqsubseteq$1 = "⊑";
var sqsup$1 = "⊐";
var sqsupe$1 = "⊒";
var sqsupset$1 = "⊐";
var sqsupseteq$1 = "⊒";
var square$1 = "□";
var Square$1 = "□";
var SquareIntersection$1 = "⊓";
var SquareSubset$1 = "⊏";
var SquareSubsetEqual$1 = "⊑";
var SquareSuperset$1 = "⊐";
var SquareSupersetEqual$1 = "⊒";
var SquareUnion$1 = "⊔";
var squarf$1 = "▪";
var squ$1 = "□";
var squf$1 = "▪";
var srarr$1 = "→";
var Sscr$1 = "𝒮";
var sscr$1 = "𝓈";
var ssetmn$1 = "∖";
var ssmile$1 = "⌣";
var sstarf$1 = "⋆";
var Star$1 = "⋆";
var star$1 = "☆";
var starf$1 = "★";
var straightepsilon$1 = "ϵ";
var straightphi$1 = "ϕ";
var strns$1 = "¯";
var sub$1 = "⊂";
var Sub$1 = "⋐";
var subdot$1 = "⪽";
var subE$1 = "⫅";
var sube$1 = "⊆";
var subedot$1 = "⫃";
var submult$1 = "⫁";
var subnE$1 = "⫋";
var subne$1 = "⊊";
var subplus$1 = "⪿";
var subrarr$1 = "⥹";
var subset$1 = "⊂";
var Subset$1 = "⋐";
var subseteq$1 = "⊆";
var subseteqq$1 = "⫅";
var SubsetEqual$1 = "⊆";
var subsetneq$1 = "⊊";
var subsetneqq$1 = "⫋";
var subsim$1 = "⫇";
var subsub$1 = "⫕";
var subsup$1 = "⫓";
var succapprox$1 = "⪸";
var succ$1 = "≻";
var succcurlyeq$1 = "≽";
var Succeeds$1 = "≻";
var SucceedsEqual$1 = "⪰";
var SucceedsSlantEqual$1 = "≽";
var SucceedsTilde$1 = "≿";
var succeq$1 = "⪰";
var succnapprox$1 = "⪺";
var succneqq$1 = "⪶";
var succnsim$1 = "⋩";
var succsim$1 = "≿";
var SuchThat$1 = "∋";
var sum$1 = "∑";
var Sum$1 = "∑";
var sung$1 = "♪";
var sup1$3 = "¹";
var sup2$3 = "²";
var sup3$3 = "³";
var sup$1 = "⊃";
var Sup$1 = "⋑";
var supdot$1 = "⪾";
var supdsub$1 = "⫘";
var supE$1 = "⫆";
var supe$1 = "⊇";
var supedot$1 = "⫄";
var Superset$1 = "⊃";
var SupersetEqual$1 = "⊇";
var suphsol$1 = "⟉";
var suphsub$1 = "⫗";
var suplarr$1 = "⥻";
var supmult$1 = "⫂";
var supnE$1 = "⫌";
var supne$1 = "⊋";
var supplus$1 = "⫀";
var supset$1 = "⊃";
var Supset$1 = "⋑";
var supseteq$1 = "⊇";
var supseteqq$1 = "⫆";
var supsetneq$1 = "⊋";
var supsetneqq$1 = "⫌";
var supsim$1 = "⫈";
var supsub$1 = "⫔";
var supsup$1 = "⫖";
var swarhk$1 = "⤦";
var swarr$1 = "↙";
var swArr$1 = "⇙";
var swarrow$1 = "↙";
var swnwar$1 = "⤪";
var szlig$3 = "ß";
var Tab$1 = "\t";
var target$1 = "⌖";
var Tau$1 = "Τ";
var tau$1 = "τ";
var tbrk$1 = "⎴";
var Tcaron$1 = "Ť";
var tcaron$1 = "ť";
var Tcedil$1 = "Ţ";
var tcedil$1 = "ţ";
var Tcy$1 = "Т";
var tcy$1 = "т";
var tdot$1 = "⃛";
var telrec$1 = "⌕";
var Tfr$1 = "𝔗";
var tfr$1 = "𝔱";
var there4$1 = "∴";
var therefore$1 = "∴";
var Therefore$1 = "∴";
var Theta$1 = "Θ";
var theta$1 = "θ";
var thetasym$1 = "ϑ";
var thetav$1 = "ϑ";
var thickapprox$1 = "≈";
var thicksim$1 = "∼";
var ThickSpace$1 = "  ";
var ThinSpace$1 = " ";
var thinsp$1 = " ";
var thkap$1 = "≈";
var thksim$1 = "∼";
var THORN$3 = "Þ";
var thorn$3 = "þ";
var tilde$2 = "˜";
var Tilde$1 = "∼";
var TildeEqual$1 = "≃";
var TildeFullEqual$1 = "≅";
var TildeTilde$1 = "≈";
var timesbar$1 = "⨱";
var timesb$1 = "⊠";
var times$3 = "×";
var timesd$1 = "⨰";
var tint$1 = "∭";
var toea$1 = "⤨";
var topbot$1 = "⌶";
var topcir$1 = "⫱";
var top$1 = "⊤";
var Topf$1 = "𝕋";
var topf$1 = "𝕥";
var topfork$1 = "⫚";
var tosa$1 = "⤩";
var tprime$1 = "‴";
var trade$1 = "™";
var TRADE$1 = "™";
var triangle$1 = "▵";
var triangledown$1 = "▿";
var triangleleft$1 = "◃";
var trianglelefteq$1 = "⊴";
var triangleq$1 = "≜";
var triangleright$1 = "▹";
var trianglerighteq$1 = "⊵";
var tridot$1 = "◬";
var trie$1 = "≜";
var triminus$1 = "⨺";
var TripleDot$1 = "⃛";
var triplus$1 = "⨹";
var trisb$1 = "⧍";
var tritime$1 = "⨻";
var trpezium$1 = "⏢";
var Tscr$1 = "𝒯";
var tscr$1 = "𝓉";
var TScy$1 = "Ц";
var tscy$1 = "ц";
var TSHcy$1 = "Ћ";
var tshcy$1 = "ћ";
var Tstrok$1 = "Ŧ";
var tstrok$1 = "ŧ";
var twixt$1 = "≬";
var twoheadleftarrow$1 = "↞";
var twoheadrightarrow$1 = "↠";
var Uacute$3 = "Ú";
var uacute$3 = "ú";
var uarr$1 = "↑";
var Uarr$1 = "↟";
var uArr$1 = "⇑";
var Uarrocir$1 = "⥉";
var Ubrcy$1 = "Ў";
var ubrcy$1 = "ў";
var Ubreve$1 = "Ŭ";
var ubreve$1 = "ŭ";
var Ucirc$3 = "Û";
var ucirc$3 = "û";
var Ucy$1 = "У";
var ucy$1 = "у";
var udarr$1 = "⇅";
var Udblac$1 = "Ű";
var udblac$1 = "ű";
var udhar$1 = "⥮";
var ufisht$1 = "⥾";
var Ufr$1 = "𝔘";
var ufr$1 = "𝔲";
var Ugrave$3 = "Ù";
var ugrave$3 = "ù";
var uHar$1 = "⥣";
var uharl$1 = "↿";
var uharr$1 = "↾";
var uhblk$1 = "▀";
var ulcorn$1 = "⌜";
var ulcorner$1 = "⌜";
var ulcrop$1 = "⌏";
var ultri$1 = "◸";
var Umacr$1 = "Ū";
var umacr$1 = "ū";
var uml$3 = "¨";
var UnderBar$1 = "_";
var UnderBrace$1 = "⏟";
var UnderBracket$1 = "⎵";
var UnderParenthesis$1 = "⏝";
var Union$1 = "⋃";
var UnionPlus$1 = "⊎";
var Uogon$1 = "Ų";
var uogon$1 = "ų";
var Uopf$1 = "𝕌";
var uopf$1 = "𝕦";
var UpArrowBar$1 = "⤒";
var uparrow$1 = "↑";
var UpArrow$1 = "↑";
var Uparrow$1 = "⇑";
var UpArrowDownArrow$1 = "⇅";
var updownarrow$1 = "↕";
var UpDownArrow$1 = "↕";
var Updownarrow$1 = "⇕";
var UpEquilibrium$1 = "⥮";
var upharpoonleft$1 = "↿";
var upharpoonright$1 = "↾";
var uplus$1 = "⊎";
var UpperLeftArrow$1 = "↖";
var UpperRightArrow$1 = "↗";
var upsi$1 = "υ";
var Upsi$1 = "ϒ";
var upsih$1 = "ϒ";
var Upsilon$1 = "Υ";
var upsilon$1 = "υ";
var UpTeeArrow$1 = "↥";
var UpTee$1 = "⊥";
var upuparrows$1 = "⇈";
var urcorn$1 = "⌝";
var urcorner$1 = "⌝";
var urcrop$1 = "⌎";
var Uring$1 = "Ů";
var uring$1 = "ů";
var urtri$1 = "◹";
var Uscr$1 = "𝒰";
var uscr$1 = "𝓊";
var utdot$1 = "⋰";
var Utilde$1 = "Ũ";
var utilde$1 = "ũ";
var utri$1 = "▵";
var utrif$1 = "▴";
var uuarr$1 = "⇈";
var Uuml$3 = "Ü";
var uuml$3 = "ü";
var uwangle$1 = "⦧";
var vangrt$1 = "⦜";
var varepsilon$1 = "ϵ";
var varkappa$1 = "ϰ";
var varnothing$1 = "∅";
var varphi$1 = "ϕ";
var varpi$1 = "ϖ";
var varpropto$1 = "∝";
var varr$1 = "↕";
var vArr$1 = "⇕";
var varrho$1 = "ϱ";
var varsigma$1 = "ς";
var varsubsetneq$1 = "⊊︀";
var varsubsetneqq$1 = "⫋︀";
var varsupsetneq$1 = "⊋︀";
var varsupsetneqq$1 = "⫌︀";
var vartheta$1 = "ϑ";
var vartriangleleft$1 = "⊲";
var vartriangleright$1 = "⊳";
var vBar$1 = "⫨";
var Vbar$1 = "⫫";
var vBarv$1 = "⫩";
var Vcy$1 = "В";
var vcy$1 = "в";
var vdash$1 = "⊢";
var vDash$1 = "⊨";
var Vdash$1 = "⊩";
var VDash$1 = "⊫";
var Vdashl$1 = "⫦";
var veebar$1 = "⊻";
var vee$1 = "∨";
var Vee$1 = "⋁";
var veeeq$1 = "≚";
var vellip$1 = "⋮";
var verbar$1 = "|";
var Verbar$1 = "‖";
var vert$1 = "|";
var Vert$1 = "‖";
var VerticalBar$1 = "∣";
var VerticalLine$1 = "|";
var VerticalSeparator$1 = "❘";
var VerticalTilde$1 = "≀";
var VeryThinSpace$1 = " ";
var Vfr$1 = "𝔙";
var vfr$1 = "𝔳";
var vltri$1 = "⊲";
var vnsub$1 = "⊂⃒";
var vnsup$1 = "⊃⃒";
var Vopf$1 = "𝕍";
var vopf$1 = "𝕧";
var vprop$1 = "∝";
var vrtri$1 = "⊳";
var Vscr$1 = "𝒱";
var vscr$1 = "𝓋";
var vsubnE$1 = "⫋︀";
var vsubne$1 = "⊊︀";
var vsupnE$1 = "⫌︀";
var vsupne$1 = "⊋︀";
var Vvdash$1 = "⊪";
var vzigzag$1 = "⦚";
var Wcirc$1 = "Ŵ";
var wcirc$1 = "ŵ";
var wedbar$1 = "⩟";
var wedge$1 = "∧";
var Wedge$1 = "⋀";
var wedgeq$1 = "≙";
var weierp$1 = "℘";
var Wfr$1 = "𝔚";
var wfr$1 = "𝔴";
var Wopf$1 = "𝕎";
var wopf$1 = "𝕨";
var wp$1 = "℘";
var wr$1 = "≀";
var wreath$1 = "≀";
var Wscr$1 = "𝒲";
var wscr$1 = "𝓌";
var xcap$1 = "⋂";
var xcirc$1 = "◯";
var xcup$1 = "⋃";
var xdtri$1 = "▽";
var Xfr$1 = "𝔛";
var xfr$1 = "𝔵";
var xharr$1 = "⟷";
var xhArr$1 = "⟺";
var Xi$1 = "Ξ";
var xi$1 = "ξ";
var xlarr$1 = "⟵";
var xlArr$1 = "⟸";
var xmap$1 = "⟼";
var xnis$1 = "⋻";
var xodot$1 = "⨀";
var Xopf$1 = "𝕏";
var xopf$1 = "𝕩";
var xoplus$1 = "⨁";
var xotime$1 = "⨂";
var xrarr$1 = "⟶";
var xrArr$1 = "⟹";
var Xscr$1 = "𝒳";
var xscr$1 = "𝓍";
var xsqcup$1 = "⨆";
var xuplus$1 = "⨄";
var xutri$1 = "△";
var xvee$1 = "⋁";
var xwedge$1 = "⋀";
var Yacute$3 = "Ý";
var yacute$3 = "ý";
var YAcy$1 = "Я";
var yacy$1 = "я";
var Ycirc$1 = "Ŷ";
var ycirc$1 = "ŷ";
var Ycy$1 = "Ы";
var ycy$1 = "ы";
var yen$3 = "¥";
var Yfr$1 = "𝔜";
var yfr$1 = "𝔶";
var YIcy$1 = "Ї";
var yicy$1 = "ї";
var Yopf$1 = "𝕐";
var yopf$1 = "𝕪";
var Yscr$1 = "𝒴";
var yscr$1 = "𝓎";
var YUcy$1 = "Ю";
var yucy$1 = "ю";
var yuml$3 = "ÿ";
var Yuml$1 = "Ÿ";
var Zacute$1 = "Ź";
var zacute$1 = "ź";
var Zcaron$1 = "Ž";
var zcaron$1 = "ž";
var Zcy$1 = "З";
var zcy$1 = "з";
var Zdot$1 = "Ż";
var zdot$1 = "ż";
var zeetrf$1 = "ℨ";
var ZeroWidthSpace$1 = "​";
var Zeta$1 = "Ζ";
var zeta$1 = "ζ";
var zfr$1 = "𝔷";
var Zfr$1 = "ℨ";
var ZHcy$1 = "Ж";
var zhcy$1 = "ж";
var zigrarr$1 = "⇝";
var zopf$1 = "𝕫";
var Zopf$1 = "ℤ";
var Zscr$1 = "𝒵";
var zscr$1 = "𝓏";
var zwj$1 = "‍";
var zwnj$1 = "‌";
var require$$1$2 = {
	Aacute: Aacute$3,
	aacute: aacute$3,
	Abreve: Abreve$1,
	abreve: abreve$1,
	ac: ac$1,
	acd: acd$1,
	acE: acE$1,
	Acirc: Acirc$3,
	acirc: acirc$3,
	acute: acute$3,
	Acy: Acy$1,
	acy: acy$1,
	AElig: AElig$3,
	aelig: aelig$3,
	af: af$1,
	Afr: Afr$1,
	afr: afr$1,
	Agrave: Agrave$3,
	agrave: agrave$3,
	alefsym: alefsym$1,
	aleph: aleph$1,
	Alpha: Alpha$1,
	alpha: alpha$1,
	Amacr: Amacr$1,
	amacr: amacr$1,
	amalg: amalg$1,
	amp: amp$5,
	AMP: AMP$3,
	andand: andand$1,
	And: And$1,
	and: and$1,
	andd: andd$1,
	andslope: andslope$1,
	andv: andv$1,
	ang: ang$1,
	ange: ange$1,
	angle: angle$1,
	angmsdaa: angmsdaa$1,
	angmsdab: angmsdab$1,
	angmsdac: angmsdac$1,
	angmsdad: angmsdad$1,
	angmsdae: angmsdae$1,
	angmsdaf: angmsdaf$1,
	angmsdag: angmsdag$1,
	angmsdah: angmsdah$1,
	angmsd: angmsd$1,
	angrt: angrt$1,
	angrtvb: angrtvb$1,
	angrtvbd: angrtvbd$1,
	angsph: angsph$1,
	angst: angst$1,
	angzarr: angzarr$1,
	Aogon: Aogon$1,
	aogon: aogon$1,
	Aopf: Aopf$1,
	aopf: aopf$1,
	apacir: apacir$1,
	ap: ap$1,
	apE: apE$1,
	ape: ape$1,
	apid: apid$1,
	apos: apos$3,
	ApplyFunction: ApplyFunction$1,
	approx: approx$1,
	approxeq: approxeq$1,
	Aring: Aring$3,
	aring: aring$3,
	Ascr: Ascr$1,
	ascr: ascr$1,
	Assign: Assign$1,
	ast: ast$1,
	asymp: asymp$1,
	asympeq: asympeq$1,
	Atilde: Atilde$3,
	atilde: atilde$3,
	Auml: Auml$3,
	auml: auml$3,
	awconint: awconint$1,
	awint: awint$1,
	backcong: backcong$1,
	backepsilon: backepsilon$1,
	backprime: backprime$1,
	backsim: backsim$1,
	backsimeq: backsimeq$1,
	Backslash: Backslash$1,
	Barv: Barv$1,
	barvee: barvee$1,
	barwed: barwed$1,
	Barwed: Barwed$1,
	barwedge: barwedge$1,
	bbrk: bbrk$1,
	bbrktbrk: bbrktbrk$1,
	bcong: bcong$1,
	Bcy: Bcy$1,
	bcy: bcy$1,
	bdquo: bdquo$1,
	becaus: becaus$1,
	because: because$1,
	Because: Because$1,
	bemptyv: bemptyv$1,
	bepsi: bepsi$1,
	bernou: bernou$1,
	Bernoullis: Bernoullis$1,
	Beta: Beta$1,
	beta: beta$1,
	beth: beth$1,
	between: between$1,
	Bfr: Bfr$1,
	bfr: bfr$1,
	bigcap: bigcap$1,
	bigcirc: bigcirc$1,
	bigcup: bigcup$1,
	bigodot: bigodot$1,
	bigoplus: bigoplus$1,
	bigotimes: bigotimes$1,
	bigsqcup: bigsqcup$1,
	bigstar: bigstar$1,
	bigtriangledown: bigtriangledown$1,
	bigtriangleup: bigtriangleup$1,
	biguplus: biguplus$1,
	bigvee: bigvee$1,
	bigwedge: bigwedge$1,
	bkarow: bkarow$1,
	blacklozenge: blacklozenge$1,
	blacksquare: blacksquare$1,
	blacktriangle: blacktriangle$1,
	blacktriangledown: blacktriangledown$1,
	blacktriangleleft: blacktriangleleft$1,
	blacktriangleright: blacktriangleright$1,
	blank: blank$1,
	blk12: blk12$1,
	blk14: blk14$1,
	blk34: blk34$1,
	block: block$1,
	bne: bne$1,
	bnequiv: bnequiv$1,
	bNot: bNot$1,
	bnot: bnot$1,
	Bopf: Bopf$1,
	bopf: bopf$1,
	bot: bot$1,
	bottom: bottom$1,
	bowtie: bowtie$1,
	boxbox: boxbox$1,
	boxdl: boxdl$1,
	boxdL: boxdL$1,
	boxDl: boxDl$1,
	boxDL: boxDL$1,
	boxdr: boxdr$1,
	boxdR: boxdR$1,
	boxDr: boxDr$1,
	boxDR: boxDR$1,
	boxh: boxh$1,
	boxH: boxH$1,
	boxhd: boxhd$1,
	boxHd: boxHd$1,
	boxhD: boxhD$1,
	boxHD: boxHD$1,
	boxhu: boxhu$1,
	boxHu: boxHu$1,
	boxhU: boxhU$1,
	boxHU: boxHU$1,
	boxminus: boxminus$1,
	boxplus: boxplus$1,
	boxtimes: boxtimes$1,
	boxul: boxul$1,
	boxuL: boxuL$1,
	boxUl: boxUl$1,
	boxUL: boxUL$1,
	boxur: boxur$1,
	boxuR: boxuR$1,
	boxUr: boxUr$1,
	boxUR: boxUR$1,
	boxv: boxv$1,
	boxV: boxV$1,
	boxvh: boxvh$1,
	boxvH: boxvH$1,
	boxVh: boxVh$1,
	boxVH: boxVH$1,
	boxvl: boxvl$1,
	boxvL: boxvL$1,
	boxVl: boxVl$1,
	boxVL: boxVL$1,
	boxvr: boxvr$1,
	boxvR: boxvR$1,
	boxVr: boxVr$1,
	boxVR: boxVR$1,
	bprime: bprime$1,
	breve: breve$1,
	Breve: Breve$1,
	brvbar: brvbar$3,
	bscr: bscr$1,
	Bscr: Bscr$1,
	bsemi: bsemi$1,
	bsim: bsim$1,
	bsime: bsime$1,
	bsolb: bsolb$1,
	bsol: bsol$1,
	bsolhsub: bsolhsub$1,
	bull: bull$1,
	bullet: bullet$1,
	bump: bump$1,
	bumpE: bumpE$1,
	bumpe: bumpe$1,
	Bumpeq: Bumpeq$1,
	bumpeq: bumpeq$1,
	Cacute: Cacute$1,
	cacute: cacute$1,
	capand: capand$1,
	capbrcup: capbrcup$1,
	capcap: capcap$1,
	cap: cap$1,
	Cap: Cap$1,
	capcup: capcup$1,
	capdot: capdot$1,
	CapitalDifferentialD: CapitalDifferentialD$1,
	caps: caps$1,
	caret: caret$2,
	caron: caron$1,
	Cayleys: Cayleys$1,
	ccaps: ccaps$1,
	Ccaron: Ccaron$1,
	ccaron: ccaron$1,
	Ccedil: Ccedil$3,
	ccedil: ccedil$3,
	Ccirc: Ccirc$1,
	ccirc: ccirc$1,
	Cconint: Cconint$1,
	ccups: ccups$1,
	ccupssm: ccupssm$1,
	Cdot: Cdot$1,
	cdot: cdot$1,
	cedil: cedil$3,
	Cedilla: Cedilla$1,
	cemptyv: cemptyv$1,
	cent: cent$3,
	centerdot: centerdot$1,
	CenterDot: CenterDot$1,
	cfr: cfr$1,
	Cfr: Cfr$1,
	CHcy: CHcy$1,
	chcy: chcy$1,
	check: check$1,
	checkmark: checkmark$1,
	Chi: Chi$1,
	chi: chi$1,
	circ: circ$1,
	circeq: circeq$1,
	circlearrowleft: circlearrowleft$1,
	circlearrowright: circlearrowright$1,
	circledast: circledast$1,
	circledcirc: circledcirc$1,
	circleddash: circleddash$1,
	CircleDot: CircleDot$1,
	circledR: circledR$1,
	circledS: circledS$1,
	CircleMinus: CircleMinus$1,
	CirclePlus: CirclePlus$1,
	CircleTimes: CircleTimes$1,
	cir: cir$1,
	cirE: cirE$1,
	cire: cire$1,
	cirfnint: cirfnint$1,
	cirmid: cirmid$1,
	cirscir: cirscir$1,
	ClockwiseContourIntegral: ClockwiseContourIntegral$1,
	CloseCurlyDoubleQuote: CloseCurlyDoubleQuote$1,
	CloseCurlyQuote: CloseCurlyQuote$1,
	clubs: clubs$1,
	clubsuit: clubsuit$1,
	colon: colon$1,
	Colon: Colon$1,
	Colone: Colone$1,
	colone: colone$1,
	coloneq: coloneq$1,
	comma: comma$2,
	commat: commat$1,
	comp: comp$1,
	compfn: compfn$1,
	complement: complement$1,
	complexes: complexes$1,
	cong: cong$1,
	congdot: congdot$1,
	Congruent: Congruent$1,
	conint: conint$1,
	Conint: Conint$1,
	ContourIntegral: ContourIntegral$1,
	copf: copf$1,
	Copf: Copf$1,
	coprod: coprod$1,
	Coproduct: Coproduct$1,
	copy: copy$3,
	COPY: COPY$3,
	copysr: copysr$1,
	CounterClockwiseContourIntegral: CounterClockwiseContourIntegral$1,
	crarr: crarr$1,
	cross: cross$1,
	Cross: Cross$1,
	Cscr: Cscr$1,
	cscr: cscr$1,
	csub: csub$1,
	csube: csube$1,
	csup: csup$1,
	csupe: csupe$1,
	ctdot: ctdot$1,
	cudarrl: cudarrl$1,
	cudarrr: cudarrr$1,
	cuepr: cuepr$1,
	cuesc: cuesc$1,
	cularr: cularr$1,
	cularrp: cularrp$1,
	cupbrcap: cupbrcap$1,
	cupcap: cupcap$1,
	CupCap: CupCap$1,
	cup: cup$1,
	Cup: Cup$1,
	cupcup: cupcup$1,
	cupdot: cupdot$1,
	cupor: cupor$1,
	cups: cups$1,
	curarr: curarr$1,
	curarrm: curarrm$1,
	curlyeqprec: curlyeqprec$1,
	curlyeqsucc: curlyeqsucc$1,
	curlyvee: curlyvee$1,
	curlywedge: curlywedge$1,
	curren: curren$3,
	curvearrowleft: curvearrowleft$1,
	curvearrowright: curvearrowright$1,
	cuvee: cuvee$1,
	cuwed: cuwed$1,
	cwconint: cwconint$1,
	cwint: cwint$1,
	cylcty: cylcty$1,
	dagger: dagger$1,
	Dagger: Dagger$1,
	daleth: daleth$1,
	darr: darr$1,
	Darr: Darr$1,
	dArr: dArr$1,
	dash: dash$1,
	Dashv: Dashv$1,
	dashv: dashv$1,
	dbkarow: dbkarow$1,
	dblac: dblac$1,
	Dcaron: Dcaron$1,
	dcaron: dcaron$1,
	Dcy: Dcy$1,
	dcy: dcy$1,
	ddagger: ddagger$1,
	ddarr: ddarr$1,
	DD: DD$1,
	dd: dd$1,
	DDotrahd: DDotrahd$1,
	ddotseq: ddotseq$1,
	deg: deg$3,
	Del: Del$1,
	Delta: Delta$1,
	delta: delta$1,
	demptyv: demptyv$1,
	dfisht: dfisht$1,
	Dfr: Dfr$1,
	dfr: dfr$1,
	dHar: dHar$1,
	dharl: dharl$1,
	dharr: dharr$1,
	DiacriticalAcute: DiacriticalAcute$1,
	DiacriticalDot: DiacriticalDot$1,
	DiacriticalDoubleAcute: DiacriticalDoubleAcute$1,
	DiacriticalGrave: DiacriticalGrave$1,
	DiacriticalTilde: DiacriticalTilde$1,
	diam: diam$1,
	diamond: diamond$1,
	Diamond: Diamond$1,
	diamondsuit: diamondsuit$1,
	diams: diams$1,
	die: die$1,
	DifferentialD: DifferentialD$1,
	digamma: digamma$1,
	disin: disin$1,
	div: div$1,
	divide: divide$3,
	divideontimes: divideontimes$1,
	divonx: divonx$1,
	DJcy: DJcy$1,
	djcy: djcy$1,
	dlcorn: dlcorn$1,
	dlcrop: dlcrop$1,
	dollar: dollar$2,
	Dopf: Dopf$1,
	dopf: dopf$1,
	Dot: Dot$1,
	dot: dot$1,
	DotDot: DotDot$1,
	doteq: doteq$1,
	doteqdot: doteqdot$1,
	DotEqual: DotEqual$1,
	dotminus: dotminus$1,
	dotplus: dotplus$1,
	dotsquare: dotsquare$1,
	doublebarwedge: doublebarwedge$1,
	DoubleContourIntegral: DoubleContourIntegral$1,
	DoubleDot: DoubleDot$1,
	DoubleDownArrow: DoubleDownArrow$1,
	DoubleLeftArrow: DoubleLeftArrow$1,
	DoubleLeftRightArrow: DoubleLeftRightArrow$1,
	DoubleLeftTee: DoubleLeftTee$1,
	DoubleLongLeftArrow: DoubleLongLeftArrow$1,
	DoubleLongLeftRightArrow: DoubleLongLeftRightArrow$1,
	DoubleLongRightArrow: DoubleLongRightArrow$1,
	DoubleRightArrow: DoubleRightArrow$1,
	DoubleRightTee: DoubleRightTee$1,
	DoubleUpArrow: DoubleUpArrow$1,
	DoubleUpDownArrow: DoubleUpDownArrow$1,
	DoubleVerticalBar: DoubleVerticalBar$1,
	DownArrowBar: DownArrowBar$1,
	downarrow: downarrow$1,
	DownArrow: DownArrow$1,
	Downarrow: Downarrow$1,
	DownArrowUpArrow: DownArrowUpArrow$1,
	DownBreve: DownBreve$1,
	downdownarrows: downdownarrows$1,
	downharpoonleft: downharpoonleft$1,
	downharpoonright: downharpoonright$1,
	DownLeftRightVector: DownLeftRightVector$1,
	DownLeftTeeVector: DownLeftTeeVector$1,
	DownLeftVectorBar: DownLeftVectorBar$1,
	DownLeftVector: DownLeftVector$1,
	DownRightTeeVector: DownRightTeeVector$1,
	DownRightVectorBar: DownRightVectorBar$1,
	DownRightVector: DownRightVector$1,
	DownTeeArrow: DownTeeArrow$1,
	DownTee: DownTee$1,
	drbkarow: drbkarow$1,
	drcorn: drcorn$1,
	drcrop: drcrop$1,
	Dscr: Dscr$1,
	dscr: dscr$1,
	DScy: DScy$1,
	dscy: dscy$1,
	dsol: dsol$1,
	Dstrok: Dstrok$1,
	dstrok: dstrok$1,
	dtdot: dtdot$1,
	dtri: dtri$1,
	dtrif: dtrif$1,
	duarr: duarr$1,
	duhar: duhar$1,
	dwangle: dwangle$1,
	DZcy: DZcy$1,
	dzcy: dzcy$1,
	dzigrarr: dzigrarr$1,
	Eacute: Eacute$3,
	eacute: eacute$3,
	easter: easter$1,
	Ecaron: Ecaron$1,
	ecaron: ecaron$1,
	Ecirc: Ecirc$3,
	ecirc: ecirc$3,
	ecir: ecir$1,
	ecolon: ecolon$1,
	Ecy: Ecy$1,
	ecy: ecy$1,
	eDDot: eDDot$1,
	Edot: Edot$1,
	edot: edot$1,
	eDot: eDot$1,
	ee: ee$1,
	efDot: efDot$1,
	Efr: Efr$1,
	efr: efr$1,
	eg: eg$1,
	Egrave: Egrave$3,
	egrave: egrave$3,
	egs: egs$1,
	egsdot: egsdot$1,
	el: el$1,
	Element: Element$1,
	elinters: elinters$1,
	ell: ell$1,
	els: els$1,
	elsdot: elsdot$1,
	Emacr: Emacr$1,
	emacr: emacr$1,
	empty: empty$1,
	emptyset: emptyset$1,
	EmptySmallSquare: EmptySmallSquare$1,
	emptyv: emptyv$1,
	EmptyVerySmallSquare: EmptyVerySmallSquare$1,
	emsp13: emsp13$1,
	emsp14: emsp14$1,
	emsp: emsp$1,
	ENG: ENG$1,
	eng: eng$1,
	ensp: ensp$1,
	Eogon: Eogon$1,
	eogon: eogon$1,
	Eopf: Eopf$1,
	eopf: eopf$1,
	epar: epar$1,
	eparsl: eparsl$1,
	eplus: eplus$1,
	epsi: epsi$1,
	Epsilon: Epsilon$1,
	epsilon: epsilon$1,
	epsiv: epsiv$1,
	eqcirc: eqcirc$1,
	eqcolon: eqcolon$1,
	eqsim: eqsim$1,
	eqslantgtr: eqslantgtr$1,
	eqslantless: eqslantless$1,
	Equal: Equal$1,
	equals: equals$1,
	EqualTilde: EqualTilde$1,
	equest: equest$1,
	Equilibrium: Equilibrium$1,
	equiv: equiv$1,
	equivDD: equivDD$1,
	eqvparsl: eqvparsl$1,
	erarr: erarr$1,
	erDot: erDot$1,
	escr: escr$1,
	Escr: Escr$1,
	esdot: esdot$1,
	Esim: Esim$1,
	esim: esim$1,
	Eta: Eta$1,
	eta: eta$1,
	ETH: ETH$3,
	eth: eth$3,
	Euml: Euml$3,
	euml: euml$3,
	euro: euro$1,
	excl: excl$1,
	exist: exist$1,
	Exists: Exists$1,
	expectation: expectation$1,
	exponentiale: exponentiale$1,
	ExponentialE: ExponentialE$1,
	fallingdotseq: fallingdotseq$1,
	Fcy: Fcy$1,
	fcy: fcy$1,
	female: female$1,
	ffilig: ffilig$1,
	fflig: fflig$1,
	ffllig: ffllig$1,
	Ffr: Ffr$1,
	ffr: ffr$1,
	filig: filig$1,
	FilledSmallSquare: FilledSmallSquare$1,
	FilledVerySmallSquare: FilledVerySmallSquare$1,
	fjlig: fjlig$1,
	flat: flat$1,
	fllig: fllig$1,
	fltns: fltns$1,
	fnof: fnof$1,
	Fopf: Fopf$1,
	fopf: fopf$1,
	forall: forall$1,
	ForAll: ForAll$1,
	fork: fork$1,
	forkv: forkv$1,
	Fouriertrf: Fouriertrf$1,
	fpartint: fpartint$1,
	frac12: frac12$3,
	frac13: frac13$1,
	frac14: frac14$3,
	frac15: frac15$1,
	frac16: frac16$1,
	frac18: frac18$1,
	frac23: frac23$1,
	frac25: frac25$1,
	frac34: frac34$3,
	frac35: frac35$1,
	frac38: frac38$1,
	frac45: frac45$1,
	frac56: frac56$1,
	frac58: frac58$1,
	frac78: frac78$1,
	frasl: frasl$1,
	frown: frown$1,
	fscr: fscr$1,
	Fscr: Fscr$1,
	gacute: gacute$1,
	Gamma: Gamma$1,
	gamma: gamma$1,
	Gammad: Gammad$1,
	gammad: gammad$1,
	gap: gap$1,
	Gbreve: Gbreve$1,
	gbreve: gbreve$1,
	Gcedil: Gcedil$1,
	Gcirc: Gcirc$1,
	gcirc: gcirc$1,
	Gcy: Gcy$1,
	gcy: gcy$1,
	Gdot: Gdot$1,
	gdot: gdot$1,
	ge: ge$1,
	gE: gE$1,
	gEl: gEl$1,
	gel: gel$1,
	geq: geq$1,
	geqq: geqq$1,
	geqslant: geqslant$1,
	gescc: gescc$1,
	ges: ges$1,
	gesdot: gesdot$1,
	gesdoto: gesdoto$1,
	gesdotol: gesdotol$1,
	gesl: gesl$1,
	gesles: gesles$1,
	Gfr: Gfr$1,
	gfr: gfr$1,
	gg: gg$1,
	Gg: Gg$1,
	ggg: ggg$1,
	gimel: gimel$1,
	GJcy: GJcy$1,
	gjcy: gjcy$1,
	gla: gla$1,
	gl: gl$1,
	glE: glE$1,
	glj: glj$1,
	gnap: gnap$1,
	gnapprox: gnapprox$1,
	gne: gne$1,
	gnE: gnE$1,
	gneq: gneq$1,
	gneqq: gneqq$1,
	gnsim: gnsim$1,
	Gopf: Gopf$1,
	gopf: gopf$1,
	grave: grave$1,
	GreaterEqual: GreaterEqual$1,
	GreaterEqualLess: GreaterEqualLess$1,
	GreaterFullEqual: GreaterFullEqual$1,
	GreaterGreater: GreaterGreater$1,
	GreaterLess: GreaterLess$1,
	GreaterSlantEqual: GreaterSlantEqual$1,
	GreaterTilde: GreaterTilde$1,
	Gscr: Gscr$1,
	gscr: gscr$1,
	gsim: gsim$1,
	gsime: gsime$1,
	gsiml: gsiml$1,
	gtcc: gtcc$1,
	gtcir: gtcir$1,
	gt: gt$6,
	GT: GT$3,
	Gt: Gt$1,
	gtdot: gtdot$1,
	gtlPar: gtlPar$1,
	gtquest: gtquest$1,
	gtrapprox: gtrapprox$1,
	gtrarr: gtrarr$1,
	gtrdot: gtrdot$1,
	gtreqless: gtreqless$1,
	gtreqqless: gtreqqless$1,
	gtrless: gtrless$1,
	gtrsim: gtrsim$1,
	gvertneqq: gvertneqq$1,
	gvnE: gvnE$1,
	Hacek: Hacek$1,
	hairsp: hairsp$1,
	half: half$1,
	hamilt: hamilt$1,
	HARDcy: HARDcy$1,
	hardcy: hardcy$1,
	harrcir: harrcir$1,
	harr: harr$1,
	hArr: hArr$1,
	harrw: harrw$1,
	Hat: Hat$1,
	hbar: hbar$1,
	Hcirc: Hcirc$1,
	hcirc: hcirc$1,
	hearts: hearts$1,
	heartsuit: heartsuit$1,
	hellip: hellip$1,
	hercon: hercon$1,
	hfr: hfr$1,
	Hfr: Hfr$1,
	HilbertSpace: HilbertSpace$1,
	hksearow: hksearow$1,
	hkswarow: hkswarow$1,
	hoarr: hoarr$1,
	homtht: homtht$1,
	hookleftarrow: hookleftarrow$1,
	hookrightarrow: hookrightarrow$1,
	hopf: hopf$1,
	Hopf: Hopf$1,
	horbar: horbar$1,
	HorizontalLine: HorizontalLine$1,
	hscr: hscr$1,
	Hscr: Hscr$1,
	hslash: hslash$1,
	Hstrok: Hstrok$1,
	hstrok: hstrok$1,
	HumpDownHump: HumpDownHump$1,
	HumpEqual: HumpEqual$1,
	hybull: hybull$1,
	hyphen: hyphen$1,
	Iacute: Iacute$3,
	iacute: iacute$3,
	ic: ic$1,
	Icirc: Icirc$3,
	icirc: icirc$3,
	Icy: Icy$1,
	icy: icy$1,
	Idot: Idot$1,
	IEcy: IEcy$1,
	iecy: iecy$1,
	iexcl: iexcl$3,
	iff: iff$1,
	ifr: ifr$1,
	Ifr: Ifr$1,
	Igrave: Igrave$3,
	igrave: igrave$3,
	ii: ii$1,
	iiiint: iiiint$1,
	iiint: iiint$1,
	iinfin: iinfin$1,
	iiota: iiota$1,
	IJlig: IJlig$1,
	ijlig: ijlig$1,
	Imacr: Imacr$1,
	imacr: imacr$1,
	image: image$1,
	ImaginaryI: ImaginaryI$1,
	imagline: imagline$1,
	imagpart: imagpart$1,
	imath: imath$1,
	Im: Im$1,
	imof: imof$1,
	imped: imped$1,
	Implies: Implies$1,
	incare: incare$1,
	"in": "∈",
	infin: infin$1,
	infintie: infintie$1,
	inodot: inodot$1,
	intcal: intcal$1,
	int: int$1,
	Int: Int$1,
	integers: integers$1,
	Integral: Integral$1,
	intercal: intercal$1,
	Intersection: Intersection$1,
	intlarhk: intlarhk$1,
	intprod: intprod$1,
	InvisibleComma: InvisibleComma$1,
	InvisibleTimes: InvisibleTimes$1,
	IOcy: IOcy$1,
	iocy: iocy$1,
	Iogon: Iogon$1,
	iogon: iogon$1,
	Iopf: Iopf$1,
	iopf: iopf$1,
	Iota: Iota$1,
	iota: iota$1,
	iprod: iprod$1,
	iquest: iquest$3,
	iscr: iscr$1,
	Iscr: Iscr$1,
	isin: isin$1,
	isindot: isindot$1,
	isinE: isinE$1,
	isins: isins$1,
	isinsv: isinsv$1,
	isinv: isinv$1,
	it: it$1,
	Itilde: Itilde$1,
	itilde: itilde$1,
	Iukcy: Iukcy$1,
	iukcy: iukcy$1,
	Iuml: Iuml$3,
	iuml: iuml$3,
	Jcirc: Jcirc$1,
	jcirc: jcirc$1,
	Jcy: Jcy$1,
	jcy: jcy$1,
	Jfr: Jfr$1,
	jfr: jfr$1,
	jmath: jmath$1,
	Jopf: Jopf$1,
	jopf: jopf$1,
	Jscr: Jscr$1,
	jscr: jscr$1,
	Jsercy: Jsercy$1,
	jsercy: jsercy$1,
	Jukcy: Jukcy$1,
	jukcy: jukcy$1,
	Kappa: Kappa$1,
	kappa: kappa$1,
	kappav: kappav$1,
	Kcedil: Kcedil$1,
	kcedil: kcedil$1,
	Kcy: Kcy$1,
	kcy: kcy$1,
	Kfr: Kfr$1,
	kfr: kfr$1,
	kgreen: kgreen$1,
	KHcy: KHcy$1,
	khcy: khcy$1,
	KJcy: KJcy$1,
	kjcy: kjcy$1,
	Kopf: Kopf$1,
	kopf: kopf$1,
	Kscr: Kscr$1,
	kscr: kscr$1,
	lAarr: lAarr$1,
	Lacute: Lacute$1,
	lacute: lacute$1,
	laemptyv: laemptyv$1,
	lagran: lagran$1,
	Lambda: Lambda$1,
	lambda: lambda$1,
	lang: lang$1,
	Lang: Lang$1,
	langd: langd$1,
	langle: langle$1,
	lap: lap$1,
	Laplacetrf: Laplacetrf$1,
	laquo: laquo$3,
	larrb: larrb$1,
	larrbfs: larrbfs$1,
	larr: larr$1,
	Larr: Larr$1,
	lArr: lArr$1,
	larrfs: larrfs$1,
	larrhk: larrhk$1,
	larrlp: larrlp$1,
	larrpl: larrpl$1,
	larrsim: larrsim$1,
	larrtl: larrtl$1,
	latail: latail$1,
	lAtail: lAtail$1,
	lat: lat$1,
	late: late$1,
	lates: lates$1,
	lbarr: lbarr$1,
	lBarr: lBarr$1,
	lbbrk: lbbrk$1,
	lbrace: lbrace$1,
	lbrack: lbrack$1,
	lbrke: lbrke$1,
	lbrksld: lbrksld$1,
	lbrkslu: lbrkslu$1,
	Lcaron: Lcaron$1,
	lcaron: lcaron$1,
	Lcedil: Lcedil$1,
	lcedil: lcedil$1,
	lceil: lceil$1,
	lcub: lcub$1,
	Lcy: Lcy$1,
	lcy: lcy$1,
	ldca: ldca$1,
	ldquo: ldquo$1,
	ldquor: ldquor$1,
	ldrdhar: ldrdhar$1,
	ldrushar: ldrushar$1,
	ldsh: ldsh$1,
	le: le$1,
	lE: lE$1,
	LeftAngleBracket: LeftAngleBracket$1,
	LeftArrowBar: LeftArrowBar$1,
	leftarrow: leftarrow$1,
	LeftArrow: LeftArrow$1,
	Leftarrow: Leftarrow$1,
	LeftArrowRightArrow: LeftArrowRightArrow$1,
	leftarrowtail: leftarrowtail$1,
	LeftCeiling: LeftCeiling$1,
	LeftDoubleBracket: LeftDoubleBracket$1,
	LeftDownTeeVector: LeftDownTeeVector$1,
	LeftDownVectorBar: LeftDownVectorBar$1,
	LeftDownVector: LeftDownVector$1,
	LeftFloor: LeftFloor$1,
	leftharpoondown: leftharpoondown$1,
	leftharpoonup: leftharpoonup$1,
	leftleftarrows: leftleftarrows$1,
	leftrightarrow: leftrightarrow$1,
	LeftRightArrow: LeftRightArrow$1,
	Leftrightarrow: Leftrightarrow$1,
	leftrightarrows: leftrightarrows$1,
	leftrightharpoons: leftrightharpoons$1,
	leftrightsquigarrow: leftrightsquigarrow$1,
	LeftRightVector: LeftRightVector$1,
	LeftTeeArrow: LeftTeeArrow$1,
	LeftTee: LeftTee$1,
	LeftTeeVector: LeftTeeVector$1,
	leftthreetimes: leftthreetimes$1,
	LeftTriangleBar: LeftTriangleBar$1,
	LeftTriangle: LeftTriangle$1,
	LeftTriangleEqual: LeftTriangleEqual$1,
	LeftUpDownVector: LeftUpDownVector$1,
	LeftUpTeeVector: LeftUpTeeVector$1,
	LeftUpVectorBar: LeftUpVectorBar$1,
	LeftUpVector: LeftUpVector$1,
	LeftVectorBar: LeftVectorBar$1,
	LeftVector: LeftVector$1,
	lEg: lEg$1,
	leg: leg$1,
	leq: leq$1,
	leqq: leqq$1,
	leqslant: leqslant$1,
	lescc: lescc$1,
	les: les$1,
	lesdot: lesdot$1,
	lesdoto: lesdoto$1,
	lesdotor: lesdotor$1,
	lesg: lesg$1,
	lesges: lesges$1,
	lessapprox: lessapprox$1,
	lessdot: lessdot$1,
	lesseqgtr: lesseqgtr$1,
	lesseqqgtr: lesseqqgtr$1,
	LessEqualGreater: LessEqualGreater$1,
	LessFullEqual: LessFullEqual$1,
	LessGreater: LessGreater$1,
	lessgtr: lessgtr$1,
	LessLess: LessLess$1,
	lesssim: lesssim$1,
	LessSlantEqual: LessSlantEqual$1,
	LessTilde: LessTilde$1,
	lfisht: lfisht$1,
	lfloor: lfloor$1,
	Lfr: Lfr$1,
	lfr: lfr$1,
	lg: lg$1,
	lgE: lgE$1,
	lHar: lHar$1,
	lhard: lhard$1,
	lharu: lharu$1,
	lharul: lharul$1,
	lhblk: lhblk$1,
	LJcy: LJcy$1,
	ljcy: ljcy$1,
	llarr: llarr$1,
	ll: ll$1,
	Ll: Ll$1,
	llcorner: llcorner$1,
	Lleftarrow: Lleftarrow$1,
	llhard: llhard$1,
	lltri: lltri$1,
	Lmidot: Lmidot$1,
	lmidot: lmidot$1,
	lmoustache: lmoustache$1,
	lmoust: lmoust$1,
	lnap: lnap$1,
	lnapprox: lnapprox$1,
	lne: lne$1,
	lnE: lnE$1,
	lneq: lneq$1,
	lneqq: lneqq$1,
	lnsim: lnsim$1,
	loang: loang$1,
	loarr: loarr$1,
	lobrk: lobrk$1,
	longleftarrow: longleftarrow$1,
	LongLeftArrow: LongLeftArrow$1,
	Longleftarrow: Longleftarrow$1,
	longleftrightarrow: longleftrightarrow$1,
	LongLeftRightArrow: LongLeftRightArrow$1,
	Longleftrightarrow: Longleftrightarrow$1,
	longmapsto: longmapsto$1,
	longrightarrow: longrightarrow$1,
	LongRightArrow: LongRightArrow$1,
	Longrightarrow: Longrightarrow$1,
	looparrowleft: looparrowleft$1,
	looparrowright: looparrowright$1,
	lopar: lopar$1,
	Lopf: Lopf$1,
	lopf: lopf$1,
	loplus: loplus$1,
	lotimes: lotimes$1,
	lowast: lowast$1,
	lowbar: lowbar$1,
	LowerLeftArrow: LowerLeftArrow$1,
	LowerRightArrow: LowerRightArrow$1,
	loz: loz$1,
	lozenge: lozenge$1,
	lozf: lozf$1,
	lpar: lpar$1,
	lparlt: lparlt$1,
	lrarr: lrarr$1,
	lrcorner: lrcorner$1,
	lrhar: lrhar$1,
	lrhard: lrhard$1,
	lrm: lrm$1,
	lrtri: lrtri$1,
	lsaquo: lsaquo$1,
	lscr: lscr$1,
	Lscr: Lscr$1,
	lsh: lsh$1,
	Lsh: Lsh$1,
	lsim: lsim$1,
	lsime: lsime$1,
	lsimg: lsimg$1,
	lsqb: lsqb$1,
	lsquo: lsquo$1,
	lsquor: lsquor$1,
	Lstrok: Lstrok$1,
	lstrok: lstrok$1,
	ltcc: ltcc$1,
	ltcir: ltcir$1,
	lt: lt$5,
	LT: LT$3,
	Lt: Lt$1,
	ltdot: ltdot$1,
	lthree: lthree$1,
	ltimes: ltimes$1,
	ltlarr: ltlarr$1,
	ltquest: ltquest$1,
	ltri: ltri$1,
	ltrie: ltrie$1,
	ltrif: ltrif$1,
	ltrPar: ltrPar$1,
	lurdshar: lurdshar$1,
	luruhar: luruhar$1,
	lvertneqq: lvertneqq$1,
	lvnE: lvnE$1,
	macr: macr$3,
	male: male$1,
	malt: malt$1,
	maltese: maltese$1,
	"Map": "⤅",
	map: map$1,
	mapsto: mapsto$1,
	mapstodown: mapstodown$1,
	mapstoleft: mapstoleft$1,
	mapstoup: mapstoup$1,
	marker: marker$1,
	mcomma: mcomma$1,
	Mcy: Mcy$1,
	mcy: mcy$1,
	mdash: mdash$1,
	mDDot: mDDot$1,
	measuredangle: measuredangle$1,
	MediumSpace: MediumSpace$1,
	Mellintrf: Mellintrf$1,
	Mfr: Mfr$1,
	mfr: mfr$1,
	mho: mho$1,
	micro: micro$3,
	midast: midast$1,
	midcir: midcir$1,
	mid: mid$1,
	middot: middot$3,
	minusb: minusb$1,
	minus: minus$1,
	minusd: minusd$1,
	minusdu: minusdu$1,
	MinusPlus: MinusPlus$1,
	mlcp: mlcp$1,
	mldr: mldr$1,
	mnplus: mnplus$1,
	models: models$1,
	Mopf: Mopf$1,
	mopf: mopf$1,
	mp: mp$1,
	mscr: mscr$1,
	Mscr: Mscr$1,
	mstpos: mstpos$1,
	Mu: Mu$1,
	mu: mu$1,
	multimap: multimap$1,
	mumap: mumap$1,
	nabla: nabla$1,
	Nacute: Nacute$1,
	nacute: nacute$1,
	nang: nang$1,
	nap: nap$1,
	napE: napE$1,
	napid: napid$1,
	napos: napos$1,
	napprox: napprox$1,
	natural: natural$1,
	naturals: naturals$1,
	natur: natur$1,
	nbsp: nbsp$3,
	nbump: nbump$1,
	nbumpe: nbumpe$1,
	ncap: ncap$1,
	Ncaron: Ncaron$1,
	ncaron: ncaron$1,
	Ncedil: Ncedil$1,
	ncedil: ncedil$1,
	ncong: ncong$1,
	ncongdot: ncongdot$1,
	ncup: ncup$1,
	Ncy: Ncy$1,
	ncy: ncy$1,
	ndash: ndash$1,
	nearhk: nearhk$1,
	nearr: nearr$1,
	neArr: neArr$1,
	nearrow: nearrow$1,
	ne: ne$1,
	nedot: nedot$1,
	NegativeMediumSpace: NegativeMediumSpace$1,
	NegativeThickSpace: NegativeThickSpace$1,
	NegativeThinSpace: NegativeThinSpace$1,
	NegativeVeryThinSpace: NegativeVeryThinSpace$1,
	nequiv: nequiv$1,
	nesear: nesear$1,
	nesim: nesim$1,
	NestedGreaterGreater: NestedGreaterGreater$1,
	NestedLessLess: NestedLessLess$1,
	NewLine: NewLine$1,
	nexist: nexist$1,
	nexists: nexists$1,
	Nfr: Nfr$1,
	nfr: nfr$1,
	ngE: ngE$1,
	nge: nge$1,
	ngeq: ngeq$1,
	ngeqq: ngeqq$1,
	ngeqslant: ngeqslant$1,
	nges: nges$1,
	nGg: nGg$1,
	ngsim: ngsim$1,
	nGt: nGt$1,
	ngt: ngt$1,
	ngtr: ngtr$1,
	nGtv: nGtv$1,
	nharr: nharr$1,
	nhArr: nhArr$1,
	nhpar: nhpar$1,
	ni: ni$1,
	nis: nis$1,
	nisd: nisd$1,
	niv: niv$1,
	NJcy: NJcy$1,
	njcy: njcy$1,
	nlarr: nlarr$1,
	nlArr: nlArr$1,
	nldr: nldr$1,
	nlE: nlE$1,
	nle: nle$1,
	nleftarrow: nleftarrow$1,
	nLeftarrow: nLeftarrow$1,
	nleftrightarrow: nleftrightarrow$1,
	nLeftrightarrow: nLeftrightarrow$1,
	nleq: nleq$1,
	nleqq: nleqq$1,
	nleqslant: nleqslant$1,
	nles: nles$1,
	nless: nless$1,
	nLl: nLl$1,
	nlsim: nlsim$1,
	nLt: nLt$1,
	nlt: nlt$1,
	nltri: nltri$1,
	nltrie: nltrie$1,
	nLtv: nLtv$1,
	nmid: nmid$1,
	NoBreak: NoBreak$1,
	NonBreakingSpace: NonBreakingSpace$1,
	nopf: nopf$1,
	Nopf: Nopf$1,
	Not: Not$1,
	not: not$3,
	NotCongruent: NotCongruent$1,
	NotCupCap: NotCupCap$1,
	NotDoubleVerticalBar: NotDoubleVerticalBar$1,
	NotElement: NotElement$1,
	NotEqual: NotEqual$1,
	NotEqualTilde: NotEqualTilde$1,
	NotExists: NotExists$1,
	NotGreater: NotGreater$1,
	NotGreaterEqual: NotGreaterEqual$1,
	NotGreaterFullEqual: NotGreaterFullEqual$1,
	NotGreaterGreater: NotGreaterGreater$1,
	NotGreaterLess: NotGreaterLess$1,
	NotGreaterSlantEqual: NotGreaterSlantEqual$1,
	NotGreaterTilde: NotGreaterTilde$1,
	NotHumpDownHump: NotHumpDownHump$1,
	NotHumpEqual: NotHumpEqual$1,
	notin: notin$1,
	notindot: notindot$1,
	notinE: notinE$1,
	notinva: notinva$1,
	notinvb: notinvb$1,
	notinvc: notinvc$1,
	NotLeftTriangleBar: NotLeftTriangleBar$1,
	NotLeftTriangle: NotLeftTriangle$1,
	NotLeftTriangleEqual: NotLeftTriangleEqual$1,
	NotLess: NotLess$1,
	NotLessEqual: NotLessEqual$1,
	NotLessGreater: NotLessGreater$1,
	NotLessLess: NotLessLess$1,
	NotLessSlantEqual: NotLessSlantEqual$1,
	NotLessTilde: NotLessTilde$1,
	NotNestedGreaterGreater: NotNestedGreaterGreater$1,
	NotNestedLessLess: NotNestedLessLess$1,
	notni: notni$1,
	notniva: notniva$1,
	notnivb: notnivb$1,
	notnivc: notnivc$1,
	NotPrecedes: NotPrecedes$1,
	NotPrecedesEqual: NotPrecedesEqual$1,
	NotPrecedesSlantEqual: NotPrecedesSlantEqual$1,
	NotReverseElement: NotReverseElement$1,
	NotRightTriangleBar: NotRightTriangleBar$1,
	NotRightTriangle: NotRightTriangle$1,
	NotRightTriangleEqual: NotRightTriangleEqual$1,
	NotSquareSubset: NotSquareSubset$1,
	NotSquareSubsetEqual: NotSquareSubsetEqual$1,
	NotSquareSuperset: NotSquareSuperset$1,
	NotSquareSupersetEqual: NotSquareSupersetEqual$1,
	NotSubset: NotSubset$1,
	NotSubsetEqual: NotSubsetEqual$1,
	NotSucceeds: NotSucceeds$1,
	NotSucceedsEqual: NotSucceedsEqual$1,
	NotSucceedsSlantEqual: NotSucceedsSlantEqual$1,
	NotSucceedsTilde: NotSucceedsTilde$1,
	NotSuperset: NotSuperset$1,
	NotSupersetEqual: NotSupersetEqual$1,
	NotTilde: NotTilde$1,
	NotTildeEqual: NotTildeEqual$1,
	NotTildeFullEqual: NotTildeFullEqual$1,
	NotTildeTilde: NotTildeTilde$1,
	NotVerticalBar: NotVerticalBar$1,
	nparallel: nparallel$1,
	npar: npar$1,
	nparsl: nparsl$1,
	npart: npart$1,
	npolint: npolint$1,
	npr: npr$1,
	nprcue: nprcue$1,
	nprec: nprec$1,
	npreceq: npreceq$1,
	npre: npre$1,
	nrarrc: nrarrc$1,
	nrarr: nrarr$1,
	nrArr: nrArr$1,
	nrarrw: nrarrw$1,
	nrightarrow: nrightarrow$1,
	nRightarrow: nRightarrow$1,
	nrtri: nrtri$1,
	nrtrie: nrtrie$1,
	nsc: nsc$1,
	nsccue: nsccue$1,
	nsce: nsce$1,
	Nscr: Nscr$1,
	nscr: nscr$1,
	nshortmid: nshortmid$1,
	nshortparallel: nshortparallel$1,
	nsim: nsim$1,
	nsime: nsime$1,
	nsimeq: nsimeq$1,
	nsmid: nsmid$1,
	nspar: nspar$1,
	nsqsube: nsqsube$1,
	nsqsupe: nsqsupe$1,
	nsub: nsub$1,
	nsubE: nsubE$1,
	nsube: nsube$1,
	nsubset: nsubset$1,
	nsubseteq: nsubseteq$1,
	nsubseteqq: nsubseteqq$1,
	nsucc: nsucc$1,
	nsucceq: nsucceq$1,
	nsup: nsup$1,
	nsupE: nsupE$1,
	nsupe: nsupe$1,
	nsupset: nsupset$1,
	nsupseteq: nsupseteq$1,
	nsupseteqq: nsupseteqq$1,
	ntgl: ntgl$1,
	Ntilde: Ntilde$3,
	ntilde: ntilde$3,
	ntlg: ntlg$1,
	ntriangleleft: ntriangleleft$1,
	ntrianglelefteq: ntrianglelefteq$1,
	ntriangleright: ntriangleright$1,
	ntrianglerighteq: ntrianglerighteq$1,
	Nu: Nu$1,
	nu: nu$1,
	num: num$1,
	numero: numero$1,
	numsp: numsp$1,
	nvap: nvap$1,
	nvdash: nvdash$1,
	nvDash: nvDash$1,
	nVdash: nVdash$1,
	nVDash: nVDash$1,
	nvge: nvge$1,
	nvgt: nvgt$1,
	nvHarr: nvHarr$1,
	nvinfin: nvinfin$1,
	nvlArr: nvlArr$1,
	nvle: nvle$1,
	nvlt: nvlt$1,
	nvltrie: nvltrie$1,
	nvrArr: nvrArr$1,
	nvrtrie: nvrtrie$1,
	nvsim: nvsim$1,
	nwarhk: nwarhk$1,
	nwarr: nwarr$1,
	nwArr: nwArr$1,
	nwarrow: nwarrow$1,
	nwnear: nwnear$1,
	Oacute: Oacute$3,
	oacute: oacute$3,
	oast: oast$1,
	Ocirc: Ocirc$3,
	ocirc: ocirc$3,
	ocir: ocir$1,
	Ocy: Ocy$1,
	ocy: ocy$1,
	odash: odash$1,
	Odblac: Odblac$1,
	odblac: odblac$1,
	odiv: odiv$1,
	odot: odot$1,
	odsold: odsold$1,
	OElig: OElig$1,
	oelig: oelig$1,
	ofcir: ofcir$1,
	Ofr: Ofr$1,
	ofr: ofr$1,
	ogon: ogon$1,
	Ograve: Ograve$3,
	ograve: ograve$3,
	ogt: ogt$1,
	ohbar: ohbar$1,
	ohm: ohm$1,
	oint: oint$1,
	olarr: olarr$1,
	olcir: olcir$1,
	olcross: olcross$1,
	oline: oline$1,
	olt: olt$1,
	Omacr: Omacr$1,
	omacr: omacr$1,
	Omega: Omega$1,
	omega: omega$1,
	Omicron: Omicron$1,
	omicron: omicron$1,
	omid: omid$1,
	ominus: ominus$1,
	Oopf: Oopf$1,
	oopf: oopf$1,
	opar: opar$1,
	OpenCurlyDoubleQuote: OpenCurlyDoubleQuote$1,
	OpenCurlyQuote: OpenCurlyQuote$1,
	operp: operp$1,
	oplus: oplus$1,
	orarr: orarr$1,
	Or: Or$1,
	or: or$1,
	ord: ord$1,
	order: order$1,
	orderof: orderof$1,
	ordf: ordf$3,
	ordm: ordm$3,
	origof: origof$1,
	oror: oror$1,
	orslope: orslope$1,
	orv: orv$1,
	oS: oS$1,
	Oscr: Oscr$1,
	oscr: oscr$1,
	Oslash: Oslash$3,
	oslash: oslash$3,
	osol: osol$1,
	Otilde: Otilde$3,
	otilde: otilde$3,
	otimesas: otimesas$1,
	Otimes: Otimes$1,
	otimes: otimes$1,
	Ouml: Ouml$3,
	ouml: ouml$3,
	ovbar: ovbar$1,
	OverBar: OverBar$1,
	OverBrace: OverBrace$1,
	OverBracket: OverBracket$1,
	OverParenthesis: OverParenthesis$1,
	para: para$3,
	parallel: parallel$1,
	par: par$1,
	parsim: parsim$1,
	parsl: parsl$1,
	part: part$1,
	PartialD: PartialD$1,
	Pcy: Pcy$1,
	pcy: pcy$1,
	percnt: percnt$1,
	period: period$1,
	permil: permil$1,
	perp: perp$1,
	pertenk: pertenk$1,
	Pfr: Pfr$1,
	pfr: pfr$1,
	Phi: Phi$1,
	phi: phi$1,
	phiv: phiv$1,
	phmmat: phmmat$1,
	phone: phone$1,
	Pi: Pi$1,
	pi: pi$1,
	pitchfork: pitchfork$1,
	piv: piv$1,
	planck: planck$1,
	planckh: planckh$1,
	plankv: plankv$1,
	plusacir: plusacir$1,
	plusb: plusb$1,
	pluscir: pluscir$1,
	plus: plus$2,
	plusdo: plusdo$1,
	plusdu: plusdu$1,
	pluse: pluse$1,
	PlusMinus: PlusMinus$1,
	plusmn: plusmn$3,
	plussim: plussim$1,
	plustwo: plustwo$1,
	pm: pm$1,
	Poincareplane: Poincareplane$1,
	pointint: pointint$1,
	popf: popf$1,
	Popf: Popf$1,
	pound: pound$3,
	prap: prap$1,
	Pr: Pr$1,
	pr: pr$1,
	prcue: prcue$1,
	precapprox: precapprox$1,
	prec: prec$1,
	preccurlyeq: preccurlyeq$1,
	Precedes: Precedes$1,
	PrecedesEqual: PrecedesEqual$1,
	PrecedesSlantEqual: PrecedesSlantEqual$1,
	PrecedesTilde: PrecedesTilde$1,
	preceq: preceq$1,
	precnapprox: precnapprox$1,
	precneqq: precneqq$1,
	precnsim: precnsim$1,
	pre: pre$1,
	prE: prE$1,
	precsim: precsim$1,
	prime: prime$1,
	Prime: Prime$1,
	primes: primes$1,
	prnap: prnap$1,
	prnE: prnE$1,
	prnsim: prnsim$1,
	prod: prod$1,
	Product: Product$1,
	profalar: profalar$1,
	profline: profline$1,
	profsurf: profsurf$1,
	prop: prop$1,
	Proportional: Proportional$1,
	Proportion: Proportion$1,
	propto: propto$1,
	prsim: prsim$1,
	prurel: prurel$1,
	Pscr: Pscr$1,
	pscr: pscr$1,
	Psi: Psi$1,
	psi: psi$1,
	puncsp: puncsp$1,
	Qfr: Qfr$1,
	qfr: qfr$1,
	qint: qint$1,
	qopf: qopf$1,
	Qopf: Qopf$1,
	qprime: qprime$1,
	Qscr: Qscr$1,
	qscr: qscr$1,
	quaternions: quaternions$1,
	quatint: quatint$1,
	quest: quest$1,
	questeq: questeq$1,
	quot: quot$5,
	QUOT: QUOT$3,
	rAarr: rAarr$1,
	race: race$1,
	Racute: Racute$1,
	racute: racute$1,
	radic: radic$1,
	raemptyv: raemptyv$1,
	rang: rang$1,
	Rang: Rang$1,
	rangd: rangd$1,
	range: range$1,
	rangle: rangle$1,
	raquo: raquo$3,
	rarrap: rarrap$1,
	rarrb: rarrb$1,
	rarrbfs: rarrbfs$1,
	rarrc: rarrc$1,
	rarr: rarr$1,
	Rarr: Rarr$1,
	rArr: rArr$1,
	rarrfs: rarrfs$1,
	rarrhk: rarrhk$1,
	rarrlp: rarrlp$1,
	rarrpl: rarrpl$1,
	rarrsim: rarrsim$1,
	Rarrtl: Rarrtl$1,
	rarrtl: rarrtl$1,
	rarrw: rarrw$1,
	ratail: ratail$1,
	rAtail: rAtail$1,
	ratio: ratio$1,
	rationals: rationals$1,
	rbarr: rbarr$1,
	rBarr: rBarr$1,
	RBarr: RBarr$1,
	rbbrk: rbbrk$1,
	rbrace: rbrace$1,
	rbrack: rbrack$1,
	rbrke: rbrke$1,
	rbrksld: rbrksld$1,
	rbrkslu: rbrkslu$1,
	Rcaron: Rcaron$1,
	rcaron: rcaron$1,
	Rcedil: Rcedil$1,
	rcedil: rcedil$1,
	rceil: rceil$1,
	rcub: rcub$1,
	Rcy: Rcy$1,
	rcy: rcy$1,
	rdca: rdca$1,
	rdldhar: rdldhar$1,
	rdquo: rdquo$1,
	rdquor: rdquor$1,
	rdsh: rdsh$1,
	real: real$1,
	realine: realine$1,
	realpart: realpart$1,
	reals: reals$1,
	Re: Re$1,
	rect: rect$1,
	reg: reg$3,
	REG: REG$3,
	ReverseElement: ReverseElement$1,
	ReverseEquilibrium: ReverseEquilibrium$1,
	ReverseUpEquilibrium: ReverseUpEquilibrium$1,
	rfisht: rfisht$1,
	rfloor: rfloor$1,
	rfr: rfr$1,
	Rfr: Rfr$1,
	rHar: rHar$1,
	rhard: rhard$1,
	rharu: rharu$1,
	rharul: rharul$1,
	Rho: Rho$1,
	rho: rho$1,
	rhov: rhov$1,
	RightAngleBracket: RightAngleBracket$1,
	RightArrowBar: RightArrowBar$1,
	rightarrow: rightarrow$1,
	RightArrow: RightArrow$1,
	Rightarrow: Rightarrow$1,
	RightArrowLeftArrow: RightArrowLeftArrow$1,
	rightarrowtail: rightarrowtail$1,
	RightCeiling: RightCeiling$1,
	RightDoubleBracket: RightDoubleBracket$1,
	RightDownTeeVector: RightDownTeeVector$1,
	RightDownVectorBar: RightDownVectorBar$1,
	RightDownVector: RightDownVector$1,
	RightFloor: RightFloor$1,
	rightharpoondown: rightharpoondown$1,
	rightharpoonup: rightharpoonup$1,
	rightleftarrows: rightleftarrows$1,
	rightleftharpoons: rightleftharpoons$1,
	rightrightarrows: rightrightarrows$1,
	rightsquigarrow: rightsquigarrow$1,
	RightTeeArrow: RightTeeArrow$1,
	RightTee: RightTee$1,
	RightTeeVector: RightTeeVector$1,
	rightthreetimes: rightthreetimes$1,
	RightTriangleBar: RightTriangleBar$1,
	RightTriangle: RightTriangle$1,
	RightTriangleEqual: RightTriangleEqual$1,
	RightUpDownVector: RightUpDownVector$1,
	RightUpTeeVector: RightUpTeeVector$1,
	RightUpVectorBar: RightUpVectorBar$1,
	RightUpVector: RightUpVector$1,
	RightVectorBar: RightVectorBar$1,
	RightVector: RightVector$1,
	ring: ring$1,
	risingdotseq: risingdotseq$1,
	rlarr: rlarr$1,
	rlhar: rlhar$1,
	rlm: rlm$1,
	rmoustache: rmoustache$1,
	rmoust: rmoust$1,
	rnmid: rnmid$1,
	roang: roang$1,
	roarr: roarr$1,
	robrk: robrk$1,
	ropar: ropar$1,
	ropf: ropf$1,
	Ropf: Ropf$1,
	roplus: roplus$1,
	rotimes: rotimes$1,
	RoundImplies: RoundImplies$1,
	rpar: rpar$1,
	rpargt: rpargt$1,
	rppolint: rppolint$1,
	rrarr: rrarr$1,
	Rrightarrow: Rrightarrow$1,
	rsaquo: rsaquo$1,
	rscr: rscr$1,
	Rscr: Rscr$1,
	rsh: rsh$1,
	Rsh: Rsh$1,
	rsqb: rsqb$1,
	rsquo: rsquo$1,
	rsquor: rsquor$1,
	rthree: rthree$1,
	rtimes: rtimes$1,
	rtri: rtri$1,
	rtrie: rtrie$1,
	rtrif: rtrif$1,
	rtriltri: rtriltri$1,
	RuleDelayed: RuleDelayed$1,
	ruluhar: ruluhar$1,
	rx: rx$1,
	Sacute: Sacute$1,
	sacute: sacute$1,
	sbquo: sbquo$1,
	scap: scap$1,
	Scaron: Scaron$1,
	scaron: scaron$1,
	Sc: Sc$1,
	sc: sc$1,
	sccue: sccue$1,
	sce: sce$1,
	scE: scE$1,
	Scedil: Scedil$1,
	scedil: scedil$1,
	Scirc: Scirc$1,
	scirc: scirc$1,
	scnap: scnap$1,
	scnE: scnE$1,
	scnsim: scnsim$1,
	scpolint: scpolint$1,
	scsim: scsim$1,
	Scy: Scy$1,
	scy: scy$1,
	sdotb: sdotb$1,
	sdot: sdot$1,
	sdote: sdote$1,
	searhk: searhk$1,
	searr: searr$1,
	seArr: seArr$1,
	searrow: searrow$1,
	sect: sect$3,
	semi: semi$1,
	seswar: seswar$1,
	setminus: setminus$1,
	setmn: setmn$1,
	sext: sext$1,
	Sfr: Sfr$1,
	sfr: sfr$1,
	sfrown: sfrown$1,
	sharp: sharp$1,
	SHCHcy: SHCHcy$1,
	shchcy: shchcy$1,
	SHcy: SHcy$1,
	shcy: shcy$1,
	ShortDownArrow: ShortDownArrow$1,
	ShortLeftArrow: ShortLeftArrow$1,
	shortmid: shortmid$1,
	shortparallel: shortparallel$1,
	ShortRightArrow: ShortRightArrow$1,
	ShortUpArrow: ShortUpArrow$1,
	shy: shy$3,
	Sigma: Sigma$1,
	sigma: sigma$1,
	sigmaf: sigmaf$1,
	sigmav: sigmav$1,
	sim: sim$1,
	simdot: simdot$1,
	sime: sime$1,
	simeq: simeq$1,
	simg: simg$1,
	simgE: simgE$1,
	siml: siml$1,
	simlE: simlE$1,
	simne: simne$1,
	simplus: simplus$1,
	simrarr: simrarr$1,
	slarr: slarr$1,
	SmallCircle: SmallCircle$1,
	smallsetminus: smallsetminus$1,
	smashp: smashp$1,
	smeparsl: smeparsl$1,
	smid: smid$1,
	smile: smile$1,
	smt: smt$1,
	smte: smte$1,
	smtes: smtes$1,
	SOFTcy: SOFTcy$1,
	softcy: softcy$1,
	solbar: solbar$1,
	solb: solb$1,
	sol: sol$1,
	Sopf: Sopf$1,
	sopf: sopf$1,
	spades: spades$1,
	spadesuit: spadesuit$1,
	spar: spar$1,
	sqcap: sqcap$1,
	sqcaps: sqcaps$1,
	sqcup: sqcup$1,
	sqcups: sqcups$1,
	Sqrt: Sqrt$1,
	sqsub: sqsub$1,
	sqsube: sqsube$1,
	sqsubset: sqsubset$1,
	sqsubseteq: sqsubseteq$1,
	sqsup: sqsup$1,
	sqsupe: sqsupe$1,
	sqsupset: sqsupset$1,
	sqsupseteq: sqsupseteq$1,
	square: square$1,
	Square: Square$1,
	SquareIntersection: SquareIntersection$1,
	SquareSubset: SquareSubset$1,
	SquareSubsetEqual: SquareSubsetEqual$1,
	SquareSuperset: SquareSuperset$1,
	SquareSupersetEqual: SquareSupersetEqual$1,
	SquareUnion: SquareUnion$1,
	squarf: squarf$1,
	squ: squ$1,
	squf: squf$1,
	srarr: srarr$1,
	Sscr: Sscr$1,
	sscr: sscr$1,
	ssetmn: ssetmn$1,
	ssmile: ssmile$1,
	sstarf: sstarf$1,
	Star: Star$1,
	star: star$1,
	starf: starf$1,
	straightepsilon: straightepsilon$1,
	straightphi: straightphi$1,
	strns: strns$1,
	sub: sub$1,
	Sub: Sub$1,
	subdot: subdot$1,
	subE: subE$1,
	sube: sube$1,
	subedot: subedot$1,
	submult: submult$1,
	subnE: subnE$1,
	subne: subne$1,
	subplus: subplus$1,
	subrarr: subrarr$1,
	subset: subset$1,
	Subset: Subset$1,
	subseteq: subseteq$1,
	subseteqq: subseteqq$1,
	SubsetEqual: SubsetEqual$1,
	subsetneq: subsetneq$1,
	subsetneqq: subsetneqq$1,
	subsim: subsim$1,
	subsub: subsub$1,
	subsup: subsup$1,
	succapprox: succapprox$1,
	succ: succ$1,
	succcurlyeq: succcurlyeq$1,
	Succeeds: Succeeds$1,
	SucceedsEqual: SucceedsEqual$1,
	SucceedsSlantEqual: SucceedsSlantEqual$1,
	SucceedsTilde: SucceedsTilde$1,
	succeq: succeq$1,
	succnapprox: succnapprox$1,
	succneqq: succneqq$1,
	succnsim: succnsim$1,
	succsim: succsim$1,
	SuchThat: SuchThat$1,
	sum: sum$1,
	Sum: Sum$1,
	sung: sung$1,
	sup1: sup1$3,
	sup2: sup2$3,
	sup3: sup3$3,
	sup: sup$1,
	Sup: Sup$1,
	supdot: supdot$1,
	supdsub: supdsub$1,
	supE: supE$1,
	supe: supe$1,
	supedot: supedot$1,
	Superset: Superset$1,
	SupersetEqual: SupersetEqual$1,
	suphsol: suphsol$1,
	suphsub: suphsub$1,
	suplarr: suplarr$1,
	supmult: supmult$1,
	supnE: supnE$1,
	supne: supne$1,
	supplus: supplus$1,
	supset: supset$1,
	Supset: Supset$1,
	supseteq: supseteq$1,
	supseteqq: supseteqq$1,
	supsetneq: supsetneq$1,
	supsetneqq: supsetneqq$1,
	supsim: supsim$1,
	supsub: supsub$1,
	supsup: supsup$1,
	swarhk: swarhk$1,
	swarr: swarr$1,
	swArr: swArr$1,
	swarrow: swarrow$1,
	swnwar: swnwar$1,
	szlig: szlig$3,
	Tab: Tab$1,
	target: target$1,
	Tau: Tau$1,
	tau: tau$1,
	tbrk: tbrk$1,
	Tcaron: Tcaron$1,
	tcaron: tcaron$1,
	Tcedil: Tcedil$1,
	tcedil: tcedil$1,
	Tcy: Tcy$1,
	tcy: tcy$1,
	tdot: tdot$1,
	telrec: telrec$1,
	Tfr: Tfr$1,
	tfr: tfr$1,
	there4: there4$1,
	therefore: therefore$1,
	Therefore: Therefore$1,
	Theta: Theta$1,
	theta: theta$1,
	thetasym: thetasym$1,
	thetav: thetav$1,
	thickapprox: thickapprox$1,
	thicksim: thicksim$1,
	ThickSpace: ThickSpace$1,
	ThinSpace: ThinSpace$1,
	thinsp: thinsp$1,
	thkap: thkap$1,
	thksim: thksim$1,
	THORN: THORN$3,
	thorn: thorn$3,
	tilde: tilde$2,
	Tilde: Tilde$1,
	TildeEqual: TildeEqual$1,
	TildeFullEqual: TildeFullEqual$1,
	TildeTilde: TildeTilde$1,
	timesbar: timesbar$1,
	timesb: timesb$1,
	times: times$3,
	timesd: timesd$1,
	tint: tint$1,
	toea: toea$1,
	topbot: topbot$1,
	topcir: topcir$1,
	top: top$1,
	Topf: Topf$1,
	topf: topf$1,
	topfork: topfork$1,
	tosa: tosa$1,
	tprime: tprime$1,
	trade: trade$1,
	TRADE: TRADE$1,
	triangle: triangle$1,
	triangledown: triangledown$1,
	triangleleft: triangleleft$1,
	trianglelefteq: trianglelefteq$1,
	triangleq: triangleq$1,
	triangleright: triangleright$1,
	trianglerighteq: trianglerighteq$1,
	tridot: tridot$1,
	trie: trie$1,
	triminus: triminus$1,
	TripleDot: TripleDot$1,
	triplus: triplus$1,
	trisb: trisb$1,
	tritime: tritime$1,
	trpezium: trpezium$1,
	Tscr: Tscr$1,
	tscr: tscr$1,
	TScy: TScy$1,
	tscy: tscy$1,
	TSHcy: TSHcy$1,
	tshcy: tshcy$1,
	Tstrok: Tstrok$1,
	tstrok: tstrok$1,
	twixt: twixt$1,
	twoheadleftarrow: twoheadleftarrow$1,
	twoheadrightarrow: twoheadrightarrow$1,
	Uacute: Uacute$3,
	uacute: uacute$3,
	uarr: uarr$1,
	Uarr: Uarr$1,
	uArr: uArr$1,
	Uarrocir: Uarrocir$1,
	Ubrcy: Ubrcy$1,
	ubrcy: ubrcy$1,
	Ubreve: Ubreve$1,
	ubreve: ubreve$1,
	Ucirc: Ucirc$3,
	ucirc: ucirc$3,
	Ucy: Ucy$1,
	ucy: ucy$1,
	udarr: udarr$1,
	Udblac: Udblac$1,
	udblac: udblac$1,
	udhar: udhar$1,
	ufisht: ufisht$1,
	Ufr: Ufr$1,
	ufr: ufr$1,
	Ugrave: Ugrave$3,
	ugrave: ugrave$3,
	uHar: uHar$1,
	uharl: uharl$1,
	uharr: uharr$1,
	uhblk: uhblk$1,
	ulcorn: ulcorn$1,
	ulcorner: ulcorner$1,
	ulcrop: ulcrop$1,
	ultri: ultri$1,
	Umacr: Umacr$1,
	umacr: umacr$1,
	uml: uml$3,
	UnderBar: UnderBar$1,
	UnderBrace: UnderBrace$1,
	UnderBracket: UnderBracket$1,
	UnderParenthesis: UnderParenthesis$1,
	Union: Union$1,
	UnionPlus: UnionPlus$1,
	Uogon: Uogon$1,
	uogon: uogon$1,
	Uopf: Uopf$1,
	uopf: uopf$1,
	UpArrowBar: UpArrowBar$1,
	uparrow: uparrow$1,
	UpArrow: UpArrow$1,
	Uparrow: Uparrow$1,
	UpArrowDownArrow: UpArrowDownArrow$1,
	updownarrow: updownarrow$1,
	UpDownArrow: UpDownArrow$1,
	Updownarrow: Updownarrow$1,
	UpEquilibrium: UpEquilibrium$1,
	upharpoonleft: upharpoonleft$1,
	upharpoonright: upharpoonright$1,
	uplus: uplus$1,
	UpperLeftArrow: UpperLeftArrow$1,
	UpperRightArrow: UpperRightArrow$1,
	upsi: upsi$1,
	Upsi: Upsi$1,
	upsih: upsih$1,
	Upsilon: Upsilon$1,
	upsilon: upsilon$1,
	UpTeeArrow: UpTeeArrow$1,
	UpTee: UpTee$1,
	upuparrows: upuparrows$1,
	urcorn: urcorn$1,
	urcorner: urcorner$1,
	urcrop: urcrop$1,
	Uring: Uring$1,
	uring: uring$1,
	urtri: urtri$1,
	Uscr: Uscr$1,
	uscr: uscr$1,
	utdot: utdot$1,
	Utilde: Utilde$1,
	utilde: utilde$1,
	utri: utri$1,
	utrif: utrif$1,
	uuarr: uuarr$1,
	Uuml: Uuml$3,
	uuml: uuml$3,
	uwangle: uwangle$1,
	vangrt: vangrt$1,
	varepsilon: varepsilon$1,
	varkappa: varkappa$1,
	varnothing: varnothing$1,
	varphi: varphi$1,
	varpi: varpi$1,
	varpropto: varpropto$1,
	varr: varr$1,
	vArr: vArr$1,
	varrho: varrho$1,
	varsigma: varsigma$1,
	varsubsetneq: varsubsetneq$1,
	varsubsetneqq: varsubsetneqq$1,
	varsupsetneq: varsupsetneq$1,
	varsupsetneqq: varsupsetneqq$1,
	vartheta: vartheta$1,
	vartriangleleft: vartriangleleft$1,
	vartriangleright: vartriangleright$1,
	vBar: vBar$1,
	Vbar: Vbar$1,
	vBarv: vBarv$1,
	Vcy: Vcy$1,
	vcy: vcy$1,
	vdash: vdash$1,
	vDash: vDash$1,
	Vdash: Vdash$1,
	VDash: VDash$1,
	Vdashl: Vdashl$1,
	veebar: veebar$1,
	vee: vee$1,
	Vee: Vee$1,
	veeeq: veeeq$1,
	vellip: vellip$1,
	verbar: verbar$1,
	Verbar: Verbar$1,
	vert: vert$1,
	Vert: Vert$1,
	VerticalBar: VerticalBar$1,
	VerticalLine: VerticalLine$1,
	VerticalSeparator: VerticalSeparator$1,
	VerticalTilde: VerticalTilde$1,
	VeryThinSpace: VeryThinSpace$1,
	Vfr: Vfr$1,
	vfr: vfr$1,
	vltri: vltri$1,
	vnsub: vnsub$1,
	vnsup: vnsup$1,
	Vopf: Vopf$1,
	vopf: vopf$1,
	vprop: vprop$1,
	vrtri: vrtri$1,
	Vscr: Vscr$1,
	vscr: vscr$1,
	vsubnE: vsubnE$1,
	vsubne: vsubne$1,
	vsupnE: vsupnE$1,
	vsupne: vsupne$1,
	Vvdash: Vvdash$1,
	vzigzag: vzigzag$1,
	Wcirc: Wcirc$1,
	wcirc: wcirc$1,
	wedbar: wedbar$1,
	wedge: wedge$1,
	Wedge: Wedge$1,
	wedgeq: wedgeq$1,
	weierp: weierp$1,
	Wfr: Wfr$1,
	wfr: wfr$1,
	Wopf: Wopf$1,
	wopf: wopf$1,
	wp: wp$1,
	wr: wr$1,
	wreath: wreath$1,
	Wscr: Wscr$1,
	wscr: wscr$1,
	xcap: xcap$1,
	xcirc: xcirc$1,
	xcup: xcup$1,
	xdtri: xdtri$1,
	Xfr: Xfr$1,
	xfr: xfr$1,
	xharr: xharr$1,
	xhArr: xhArr$1,
	Xi: Xi$1,
	xi: xi$1,
	xlarr: xlarr$1,
	xlArr: xlArr$1,
	xmap: xmap$1,
	xnis: xnis$1,
	xodot: xodot$1,
	Xopf: Xopf$1,
	xopf: xopf$1,
	xoplus: xoplus$1,
	xotime: xotime$1,
	xrarr: xrarr$1,
	xrArr: xrArr$1,
	Xscr: Xscr$1,
	xscr: xscr$1,
	xsqcup: xsqcup$1,
	xuplus: xuplus$1,
	xutri: xutri$1,
	xvee: xvee$1,
	xwedge: xwedge$1,
	Yacute: Yacute$3,
	yacute: yacute$3,
	YAcy: YAcy$1,
	yacy: yacy$1,
	Ycirc: Ycirc$1,
	ycirc: ycirc$1,
	Ycy: Ycy$1,
	ycy: ycy$1,
	yen: yen$3,
	Yfr: Yfr$1,
	yfr: yfr$1,
	YIcy: YIcy$1,
	yicy: yicy$1,
	Yopf: Yopf$1,
	yopf: yopf$1,
	Yscr: Yscr$1,
	yscr: yscr$1,
	YUcy: YUcy$1,
	yucy: yucy$1,
	yuml: yuml$3,
	Yuml: Yuml$1,
	Zacute: Zacute$1,
	zacute: zacute$1,
	Zcaron: Zcaron$1,
	zcaron: zcaron$1,
	Zcy: Zcy$1,
	zcy: zcy$1,
	Zdot: Zdot$1,
	zdot: zdot$1,
	zeetrf: zeetrf$1,
	ZeroWidthSpace: ZeroWidthSpace$1,
	Zeta: Zeta$1,
	zeta: zeta$1,
	zfr: zfr$1,
	Zfr: Zfr$1,
	ZHcy: ZHcy$1,
	zhcy: zhcy$1,
	zigrarr: zigrarr$1,
	zopf: zopf$1,
	Zopf: Zopf$1,
	Zscr: Zscr$1,
	zscr: zscr$1,
	zwj: zwj$1,
	zwnj: zwnj$1
};

var Aacute$2 = "Á";
var aacute$2 = "á";
var Acirc$2 = "Â";
var acirc$2 = "â";
var acute$2 = "´";
var AElig$2 = "Æ";
var aelig$2 = "æ";
var Agrave$2 = "À";
var agrave$2 = "à";
var amp$4 = "&";
var AMP$2 = "&";
var Aring$2 = "Å";
var aring$2 = "å";
var Atilde$2 = "Ã";
var atilde$2 = "ã";
var Auml$2 = "Ä";
var auml$2 = "ä";
var brvbar$2 = "¦";
var Ccedil$2 = "Ç";
var ccedil$2 = "ç";
var cedil$2 = "¸";
var cent$2 = "¢";
var copy$2 = "©";
var COPY$2 = "©";
var curren$2 = "¤";
var deg$2 = "°";
var divide$2 = "÷";
var Eacute$2 = "É";
var eacute$2 = "é";
var Ecirc$2 = "Ê";
var ecirc$2 = "ê";
var Egrave$2 = "È";
var egrave$2 = "è";
var ETH$2 = "Ð";
var eth$2 = "ð";
var Euml$2 = "Ë";
var euml$2 = "ë";
var frac12$2 = "½";
var frac14$2 = "¼";
var frac34$2 = "¾";
var gt$5 = ">";
var GT$2 = ">";
var Iacute$2 = "Í";
var iacute$2 = "í";
var Icirc$2 = "Î";
var icirc$2 = "î";
var iexcl$2 = "¡";
var Igrave$2 = "Ì";
var igrave$2 = "ì";
var iquest$2 = "¿";
var Iuml$2 = "Ï";
var iuml$2 = "ï";
var laquo$2 = "«";
var lt$4 = "<";
var LT$2 = "<";
var macr$2 = "¯";
var micro$2 = "µ";
var middot$2 = "·";
var nbsp$2 = " ";
var not$2 = "¬";
var Ntilde$2 = "Ñ";
var ntilde$2 = "ñ";
var Oacute$2 = "Ó";
var oacute$2 = "ó";
var Ocirc$2 = "Ô";
var ocirc$2 = "ô";
var Ograve$2 = "Ò";
var ograve$2 = "ò";
var ordf$2 = "ª";
var ordm$2 = "º";
var Oslash$2 = "Ø";
var oslash$2 = "ø";
var Otilde$2 = "Õ";
var otilde$2 = "õ";
var Ouml$2 = "Ö";
var ouml$2 = "ö";
var para$2 = "¶";
var plusmn$2 = "±";
var pound$2 = "£";
var quot$4 = "\"";
var QUOT$2 = "\"";
var raquo$2 = "»";
var reg$2 = "®";
var REG$2 = "®";
var sect$2 = "§";
var shy$2 = "­";
var sup1$2 = "¹";
var sup2$2 = "²";
var sup3$2 = "³";
var szlig$2 = "ß";
var THORN$2 = "Þ";
var thorn$2 = "þ";
var times$2 = "×";
var Uacute$2 = "Ú";
var uacute$2 = "ú";
var Ucirc$2 = "Û";
var ucirc$2 = "û";
var Ugrave$2 = "Ù";
var ugrave$2 = "ù";
var uml$2 = "¨";
var Uuml$2 = "Ü";
var uuml$2 = "ü";
var Yacute$2 = "Ý";
var yacute$2 = "ý";
var yen$2 = "¥";
var yuml$2 = "ÿ";
var require$$2 = {
	Aacute: Aacute$2,
	aacute: aacute$2,
	Acirc: Acirc$2,
	acirc: acirc$2,
	acute: acute$2,
	AElig: AElig$2,
	aelig: aelig$2,
	Agrave: Agrave$2,
	agrave: agrave$2,
	amp: amp$4,
	AMP: AMP$2,
	Aring: Aring$2,
	aring: aring$2,
	Atilde: Atilde$2,
	atilde: atilde$2,
	Auml: Auml$2,
	auml: auml$2,
	brvbar: brvbar$2,
	Ccedil: Ccedil$2,
	ccedil: ccedil$2,
	cedil: cedil$2,
	cent: cent$2,
	copy: copy$2,
	COPY: COPY$2,
	curren: curren$2,
	deg: deg$2,
	divide: divide$2,
	Eacute: Eacute$2,
	eacute: eacute$2,
	Ecirc: Ecirc$2,
	ecirc: ecirc$2,
	Egrave: Egrave$2,
	egrave: egrave$2,
	ETH: ETH$2,
	eth: eth$2,
	Euml: Euml$2,
	euml: euml$2,
	frac12: frac12$2,
	frac14: frac14$2,
	frac34: frac34$2,
	gt: gt$5,
	GT: GT$2,
	Iacute: Iacute$2,
	iacute: iacute$2,
	Icirc: Icirc$2,
	icirc: icirc$2,
	iexcl: iexcl$2,
	Igrave: Igrave$2,
	igrave: igrave$2,
	iquest: iquest$2,
	Iuml: Iuml$2,
	iuml: iuml$2,
	laquo: laquo$2,
	lt: lt$4,
	LT: LT$2,
	macr: macr$2,
	micro: micro$2,
	middot: middot$2,
	nbsp: nbsp$2,
	not: not$2,
	Ntilde: Ntilde$2,
	ntilde: ntilde$2,
	Oacute: Oacute$2,
	oacute: oacute$2,
	Ocirc: Ocirc$2,
	ocirc: ocirc$2,
	Ograve: Ograve$2,
	ograve: ograve$2,
	ordf: ordf$2,
	ordm: ordm$2,
	Oslash: Oslash$2,
	oslash: oslash$2,
	Otilde: Otilde$2,
	otilde: otilde$2,
	Ouml: Ouml$2,
	ouml: ouml$2,
	para: para$2,
	plusmn: plusmn$2,
	pound: pound$2,
	quot: quot$4,
	QUOT: QUOT$2,
	raquo: raquo$2,
	reg: reg$2,
	REG: REG$2,
	sect: sect$2,
	shy: shy$2,
	sup1: sup1$2,
	sup2: sup2$2,
	sup3: sup3$2,
	szlig: szlig$2,
	THORN: THORN$2,
	thorn: thorn$2,
	times: times$2,
	Uacute: Uacute$2,
	uacute: uacute$2,
	Ucirc: Ucirc$2,
	ucirc: ucirc$2,
	Ugrave: Ugrave$2,
	ugrave: ugrave$2,
	uml: uml$2,
	Uuml: Uuml$2,
	uuml: uuml$2,
	Yacute: Yacute$2,
	yacute: yacute$2,
	yen: yen$2,
	yuml: yuml$2
};

var amp$3 = "&";
var apos$2 = "'";
var gt$4 = ">";
var lt$3 = "<";
var quot$3 = "\"";
var require$$3 = {
	amp: amp$3,
	apos: apos$2,
	gt: gt$4,
	lt: lt$3,
	quot: quot$3
};

var __importDefault$6 = (commonjsGlobal && commonjsGlobal.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(Tokenizer$1, "__esModule", { value: true });
var decode_codepoint_1$1 = __importDefault$6(decode_codepoint$1);
var entities_json_1$2 = __importDefault$6(require$$1$2);
var legacy_json_1$1 = __importDefault$6(require$$2);
var xml_json_1$2 = __importDefault$6(require$$3);
function whitespace(c) {
    return c === " " || c === "\n" || c === "\t" || c === "\f" || c === "\r";
}
function isASCIIAlpha(c) {
    return (c >= "a" && c <= "z") || (c >= "A" && c <= "Z");
}
function ifElseState(upper, SUCCESS, FAILURE) {
    var lower = upper.toLowerCase();
    if (upper === lower) {
        return function (t, c) {
            if (c === lower) {
                t._state = SUCCESS;
            }
            else {
                t._state = FAILURE;
                t._index--;
            }
        };
    }
    return function (t, c) {
        if (c === lower || c === upper) {
            t._state = SUCCESS;
        }
        else {
            t._state = FAILURE;
            t._index--;
        }
    };
}
function consumeSpecialNameChar(upper, NEXT_STATE) {
    var lower = upper.toLowerCase();
    return function (t, c) {
        if (c === lower || c === upper) {
            t._state = NEXT_STATE;
        }
        else {
            t._state = 3 /* InTagName */;
            t._index--; // Consume the token again
        }
    };
}
var stateBeforeCdata1 = ifElseState("C", 24 /* BeforeCdata2 */, 16 /* InDeclaration */);
var stateBeforeCdata2 = ifElseState("D", 25 /* BeforeCdata3 */, 16 /* InDeclaration */);
var stateBeforeCdata3 = ifElseState("A", 26 /* BeforeCdata4 */, 16 /* InDeclaration */);
var stateBeforeCdata4 = ifElseState("T", 27 /* BeforeCdata5 */, 16 /* InDeclaration */);
var stateBeforeCdata5 = ifElseState("A", 28 /* BeforeCdata6 */, 16 /* InDeclaration */);
var stateBeforeScript1 = consumeSpecialNameChar("R", 35 /* BeforeScript2 */);
var stateBeforeScript2 = consumeSpecialNameChar("I", 36 /* BeforeScript3 */);
var stateBeforeScript3 = consumeSpecialNameChar("P", 37 /* BeforeScript4 */);
var stateBeforeScript4 = consumeSpecialNameChar("T", 38 /* BeforeScript5 */);
var stateAfterScript1 = ifElseState("R", 40 /* AfterScript2 */, 1 /* Text */);
var stateAfterScript2 = ifElseState("I", 41 /* AfterScript3 */, 1 /* Text */);
var stateAfterScript3 = ifElseState("P", 42 /* AfterScript4 */, 1 /* Text */);
var stateAfterScript4 = ifElseState("T", 43 /* AfterScript5 */, 1 /* Text */);
var stateBeforeStyle1 = consumeSpecialNameChar("Y", 45 /* BeforeStyle2 */);
var stateBeforeStyle2 = consumeSpecialNameChar("L", 46 /* BeforeStyle3 */);
var stateBeforeStyle3 = consumeSpecialNameChar("E", 47 /* BeforeStyle4 */);
var stateAfterStyle1 = ifElseState("Y", 49 /* AfterStyle2 */, 1 /* Text */);
var stateAfterStyle2 = ifElseState("L", 50 /* AfterStyle3 */, 1 /* Text */);
var stateAfterStyle3 = ifElseState("E", 51 /* AfterStyle4 */, 1 /* Text */);
var stateBeforeSpecialT = consumeSpecialNameChar("I", 54 /* BeforeTitle1 */);
var stateBeforeTitle1 = consumeSpecialNameChar("T", 55 /* BeforeTitle2 */);
var stateBeforeTitle2 = consumeSpecialNameChar("L", 56 /* BeforeTitle3 */);
var stateBeforeTitle3 = consumeSpecialNameChar("E", 57 /* BeforeTitle4 */);
var stateAfterSpecialTEnd = ifElseState("I", 58 /* AfterTitle1 */, 1 /* Text */);
var stateAfterTitle1 = ifElseState("T", 59 /* AfterTitle2 */, 1 /* Text */);
var stateAfterTitle2 = ifElseState("L", 60 /* AfterTitle3 */, 1 /* Text */);
var stateAfterTitle3 = ifElseState("E", 61 /* AfterTitle4 */, 1 /* Text */);
var stateBeforeEntity = ifElseState("#", 63 /* BeforeNumericEntity */, 64 /* InNamedEntity */);
var stateBeforeNumericEntity = ifElseState("X", 66 /* InHexEntity */, 65 /* InNumericEntity */);
var Tokenizer = /** @class */ (function () {
    function Tokenizer(options, cbs) {
        var _a;
        /** The current state the tokenizer is in. */
        this._state = 1 /* Text */;
        /** The read buffer. */
        this.buffer = "";
        /** The beginning of the section that is currently being read. */
        this.sectionStart = 0;
        /** The index within the buffer that we are currently looking at. */
        this._index = 0;
        /**
         * Data that has already been processed will be removed from the buffer occasionally.
         * `_bufferOffset` keeps track of how many characters have been removed, to make sure position information is accurate.
         */
        this.bufferOffset = 0;
        /** Some behavior, eg. when decoding entities, is done while we are in another state. This keeps track of the other state type. */
        this.baseState = 1 /* Text */;
        /** For special parsing behavior inside of script and style tags. */
        this.special = 1 /* None */;
        /** Indicates whether the tokenizer has been paused. */
        this.running = true;
        /** Indicates whether the tokenizer has finished running / `.end` has been called. */
        this.ended = false;
        this.cbs = cbs;
        this.xmlMode = !!(options === null || options === void 0 ? void 0 : options.xmlMode);
        this.decodeEntities = (_a = options === null || options === void 0 ? void 0 : options.decodeEntities) !== null && _a !== void 0 ? _a : true;
    }
    Tokenizer.prototype.reset = function () {
        this._state = 1 /* Text */;
        this.buffer = "";
        this.sectionStart = 0;
        this._index = 0;
        this.bufferOffset = 0;
        this.baseState = 1 /* Text */;
        this.special = 1 /* None */;
        this.running = true;
        this.ended = false;
    };
    Tokenizer.prototype.write = function (chunk) {
        if (this.ended)
            this.cbs.onerror(Error(".write() after done!"));
        this.buffer += chunk;
        this.parse();
    };
    Tokenizer.prototype.end = function (chunk) {
        if (this.ended)
            this.cbs.onerror(Error(".end() after done!"));
        if (chunk)
            this.write(chunk);
        this.ended = true;
        if (this.running)
            this.finish();
    };
    Tokenizer.prototype.pause = function () {
        this.running = false;
    };
    Tokenizer.prototype.resume = function () {
        this.running = true;
        if (this._index < this.buffer.length) {
            this.parse();
        }
        if (this.ended) {
            this.finish();
        }
    };
    /**
     * The current index within all of the written data.
     */
    Tokenizer.prototype.getAbsoluteIndex = function () {
        return this.bufferOffset + this._index;
    };
    Tokenizer.prototype.stateText = function (c) {
        if (c === "<") {
            if (this._index > this.sectionStart) {
                this.cbs.ontext(this.getSection());
            }
            this._state = 2 /* BeforeTagName */;
            this.sectionStart = this._index;
        }
        else if (this.decodeEntities &&
            c === "&" &&
            (this.special === 1 /* None */ || this.special === 4 /* Title */)) {
            if (this._index > this.sectionStart) {
                this.cbs.ontext(this.getSection());
            }
            this.baseState = 1 /* Text */;
            this._state = 62 /* BeforeEntity */;
            this.sectionStart = this._index;
        }
    };
    /**
     * HTML only allows ASCII alpha characters (a-z and A-Z) at the beginning of a tag name.
     *
     * XML allows a lot more characters here (@see https://www.w3.org/TR/REC-xml/#NT-NameStartChar).
     * We allow anything that wouldn't end the tag.
     */
    Tokenizer.prototype.isTagStartChar = function (c) {
        return (isASCIIAlpha(c) ||
            (this.xmlMode && !whitespace(c) && c !== "/" && c !== ">"));
    };
    Tokenizer.prototype.stateBeforeTagName = function (c) {
        if (c === "/") {
            this._state = 5 /* BeforeClosingTagName */;
        }
        else if (c === "<") {
            this.cbs.ontext(this.getSection());
            this.sectionStart = this._index;
        }
        else if (c === ">" ||
            this.special !== 1 /* None */ ||
            whitespace(c)) {
            this._state = 1 /* Text */;
        }
        else if (c === "!") {
            this._state = 15 /* BeforeDeclaration */;
            this.sectionStart = this._index + 1;
        }
        else if (c === "?") {
            this._state = 17 /* InProcessingInstruction */;
            this.sectionStart = this._index + 1;
        }
        else if (!this.isTagStartChar(c)) {
            this._state = 1 /* Text */;
        }
        else {
            this._state =
                !this.xmlMode && (c === "s" || c === "S")
                    ? 32 /* BeforeSpecialS */
                    : !this.xmlMode && (c === "t" || c === "T")
                        ? 52 /* BeforeSpecialT */
                        : 3 /* InTagName */;
            this.sectionStart = this._index;
        }
    };
    Tokenizer.prototype.stateInTagName = function (c) {
        if (c === "/" || c === ">" || whitespace(c)) {
            this.emitToken("onopentagname");
            this._state = 8 /* BeforeAttributeName */;
            this._index--;
        }
    };
    Tokenizer.prototype.stateBeforeClosingTagName = function (c) {
        if (whitespace(c)) ;
        else if (c === ">") {
            this._state = 1 /* Text */;
        }
        else if (this.special !== 1 /* None */) {
            if (this.special !== 4 /* Title */ && (c === "s" || c === "S")) {
                this._state = 33 /* BeforeSpecialSEnd */;
            }
            else if (this.special === 4 /* Title */ &&
                (c === "t" || c === "T")) {
                this._state = 53 /* BeforeSpecialTEnd */;
            }
            else {
                this._state = 1 /* Text */;
                this._index--;
            }
        }
        else if (!this.isTagStartChar(c)) {
            this._state = 20 /* InSpecialComment */;
            this.sectionStart = this._index;
        }
        else {
            this._state = 6 /* InClosingTagName */;
            this.sectionStart = this._index;
        }
    };
    Tokenizer.prototype.stateInClosingTagName = function (c) {
        if (c === ">" || whitespace(c)) {
            this.emitToken("onclosetag");
            this._state = 7 /* AfterClosingTagName */;
            this._index--;
        }
    };
    Tokenizer.prototype.stateAfterClosingTagName = function (c) {
        // Skip everything until ">"
        if (c === ">") {
            this._state = 1 /* Text */;
            this.sectionStart = this._index + 1;
        }
    };
    Tokenizer.prototype.stateBeforeAttributeName = function (c) {
        if (c === ">") {
            this.cbs.onopentagend();
            this._state = 1 /* Text */;
            this.sectionStart = this._index + 1;
        }
        else if (c === "/") {
            this._state = 4 /* InSelfClosingTag */;
        }
        else if (!whitespace(c)) {
            this._state = 9 /* InAttributeName */;
            this.sectionStart = this._index;
        }
    };
    Tokenizer.prototype.stateInSelfClosingTag = function (c) {
        if (c === ">") {
            this.cbs.onselfclosingtag();
            this._state = 1 /* Text */;
            this.sectionStart = this._index + 1;
            this.special = 1 /* None */; // Reset special state, in case of self-closing special tags
        }
        else if (!whitespace(c)) {
            this._state = 8 /* BeforeAttributeName */;
            this._index--;
        }
    };
    Tokenizer.prototype.stateInAttributeName = function (c) {
        if (c === "=" || c === "/" || c === ">" || whitespace(c)) {
            this.cbs.onattribname(this.getSection());
            this.sectionStart = -1;
            this._state = 10 /* AfterAttributeName */;
            this._index--;
        }
    };
    Tokenizer.prototype.stateAfterAttributeName = function (c) {
        if (c === "=") {
            this._state = 11 /* BeforeAttributeValue */;
        }
        else if (c === "/" || c === ">") {
            this.cbs.onattribend(undefined);
            this._state = 8 /* BeforeAttributeName */;
            this._index--;
        }
        else if (!whitespace(c)) {
            this.cbs.onattribend(undefined);
            this._state = 9 /* InAttributeName */;
            this.sectionStart = this._index;
        }
    };
    Tokenizer.prototype.stateBeforeAttributeValue = function (c) {
        if (c === '"') {
            this._state = 12 /* InAttributeValueDq */;
            this.sectionStart = this._index + 1;
        }
        else if (c === "'") {
            this._state = 13 /* InAttributeValueSq */;
            this.sectionStart = this._index + 1;
        }
        else if (!whitespace(c)) {
            this._state = 14 /* InAttributeValueNq */;
            this.sectionStart = this._index;
            this._index--; // Reconsume token
        }
    };
    Tokenizer.prototype.handleInAttributeValue = function (c, quote) {
        if (c === quote) {
            this.emitToken("onattribdata");
            this.cbs.onattribend(quote);
            this._state = 8 /* BeforeAttributeName */;
        }
        else if (this.decodeEntities && c === "&") {
            this.emitToken("onattribdata");
            this.baseState = this._state;
            this._state = 62 /* BeforeEntity */;
            this.sectionStart = this._index;
        }
    };
    Tokenizer.prototype.stateInAttributeValueDoubleQuotes = function (c) {
        this.handleInAttributeValue(c, '"');
    };
    Tokenizer.prototype.stateInAttributeValueSingleQuotes = function (c) {
        this.handleInAttributeValue(c, "'");
    };
    Tokenizer.prototype.stateInAttributeValueNoQuotes = function (c) {
        if (whitespace(c) || c === ">") {
            this.emitToken("onattribdata");
            this.cbs.onattribend(null);
            this._state = 8 /* BeforeAttributeName */;
            this._index--;
        }
        else if (this.decodeEntities && c === "&") {
            this.emitToken("onattribdata");
            this.baseState = this._state;
            this._state = 62 /* BeforeEntity */;
            this.sectionStart = this._index;
        }
    };
    Tokenizer.prototype.stateBeforeDeclaration = function (c) {
        this._state =
            c === "["
                ? 23 /* BeforeCdata1 */
                : c === "-"
                    ? 18 /* BeforeComment */
                    : 16 /* InDeclaration */;
    };
    Tokenizer.prototype.stateInDeclaration = function (c) {
        if (c === ">") {
            this.cbs.ondeclaration(this.getSection());
            this._state = 1 /* Text */;
            this.sectionStart = this._index + 1;
        }
    };
    Tokenizer.prototype.stateInProcessingInstruction = function (c) {
        if (c === ">") {
            this.cbs.onprocessinginstruction(this.getSection());
            this._state = 1 /* Text */;
            this.sectionStart = this._index + 1;
        }
    };
    Tokenizer.prototype.stateBeforeComment = function (c) {
        if (c === "-") {
            this._state = 19 /* InComment */;
            this.sectionStart = this._index + 1;
        }
        else {
            this._state = 16 /* InDeclaration */;
        }
    };
    Tokenizer.prototype.stateInComment = function (c) {
        if (c === "-")
            this._state = 21 /* AfterComment1 */;
    };
    Tokenizer.prototype.stateInSpecialComment = function (c) {
        if (c === ">") {
            this.cbs.oncomment(this.buffer.substring(this.sectionStart, this._index));
            this._state = 1 /* Text */;
            this.sectionStart = this._index + 1;
        }
    };
    Tokenizer.prototype.stateAfterComment1 = function (c) {
        if (c === "-") {
            this._state = 22 /* AfterComment2 */;
        }
        else {
            this._state = 19 /* InComment */;
        }
    };
    Tokenizer.prototype.stateAfterComment2 = function (c) {
        if (c === ">") {
            // Remove 2 trailing chars
            this.cbs.oncomment(this.buffer.substring(this.sectionStart, this._index - 2));
            this._state = 1 /* Text */;
            this.sectionStart = this._index + 1;
        }
        else if (c !== "-") {
            this._state = 19 /* InComment */;
        }
        // Else: stay in AFTER_COMMENT_2 (`--->`)
    };
    Tokenizer.prototype.stateBeforeCdata6 = function (c) {
        if (c === "[") {
            this._state = 29 /* InCdata */;
            this.sectionStart = this._index + 1;
        }
        else {
            this._state = 16 /* InDeclaration */;
            this._index--;
        }
    };
    Tokenizer.prototype.stateInCdata = function (c) {
        if (c === "]")
            this._state = 30 /* AfterCdata1 */;
    };
    Tokenizer.prototype.stateAfterCdata1 = function (c) {
        if (c === "]")
            this._state = 31 /* AfterCdata2 */;
        else
            this._state = 29 /* InCdata */;
    };
    Tokenizer.prototype.stateAfterCdata2 = function (c) {
        if (c === ">") {
            // Remove 2 trailing chars
            this.cbs.oncdata(this.buffer.substring(this.sectionStart, this._index - 2));
            this._state = 1 /* Text */;
            this.sectionStart = this._index + 1;
        }
        else if (c !== "]") {
            this._state = 29 /* InCdata */;
        }
        // Else: stay in AFTER_CDATA_2 (`]]]>`)
    };
    Tokenizer.prototype.stateBeforeSpecialS = function (c) {
        if (c === "c" || c === "C") {
            this._state = 34 /* BeforeScript1 */;
        }
        else if (c === "t" || c === "T") {
            this._state = 44 /* BeforeStyle1 */;
        }
        else {
            this._state = 3 /* InTagName */;
            this._index--; // Consume the token again
        }
    };
    Tokenizer.prototype.stateBeforeSpecialSEnd = function (c) {
        if (this.special === 2 /* Script */ && (c === "c" || c === "C")) {
            this._state = 39 /* AfterScript1 */;
        }
        else if (this.special === 3 /* Style */ && (c === "t" || c === "T")) {
            this._state = 48 /* AfterStyle1 */;
        }
        else
            this._state = 1 /* Text */;
    };
    Tokenizer.prototype.stateBeforeSpecialLast = function (c, special) {
        if (c === "/" || c === ">" || whitespace(c)) {
            this.special = special;
        }
        this._state = 3 /* InTagName */;
        this._index--; // Consume the token again
    };
    Tokenizer.prototype.stateAfterSpecialLast = function (c, sectionStartOffset) {
        if (c === ">" || whitespace(c)) {
            this.special = 1 /* None */;
            this._state = 6 /* InClosingTagName */;
            this.sectionStart = this._index - sectionStartOffset;
            this._index--; // Reconsume the token
        }
        else
            this._state = 1 /* Text */;
    };
    // For entities terminated with a semicolon
    Tokenizer.prototype.parseFixedEntity = function (map) {
        if (map === void 0) { map = this.xmlMode ? xml_json_1$2.default : entities_json_1$2.default; }
        // Offset = 1
        if (this.sectionStart + 1 < this._index) {
            var entity = this.buffer.substring(this.sectionStart + 1, this._index);
            if (Object.prototype.hasOwnProperty.call(map, entity)) {
                this.emitPartial(map[entity]);
                this.sectionStart = this._index + 1;
            }
        }
    };
    // Parses legacy entities (without trailing semicolon)
    Tokenizer.prototype.parseLegacyEntity = function () {
        var start = this.sectionStart + 1;
        // The max length of legacy entities is 6
        var limit = Math.min(this._index - start, 6);
        while (limit >= 2) {
            // The min length of legacy entities is 2
            var entity = this.buffer.substr(start, limit);
            if (Object.prototype.hasOwnProperty.call(legacy_json_1$1.default, entity)) {
                this.emitPartial(legacy_json_1$1.default[entity]);
                this.sectionStart += limit + 1;
                return;
            }
            limit--;
        }
    };
    Tokenizer.prototype.stateInNamedEntity = function (c) {
        if (c === ";") {
            this.parseFixedEntity();
            // Retry as legacy entity if entity wasn't parsed
            if (this.baseState === 1 /* Text */ &&
                this.sectionStart + 1 < this._index &&
                !this.xmlMode) {
                this.parseLegacyEntity();
            }
            this._state = this.baseState;
        }
        else if ((c < "0" || c > "9") && !isASCIIAlpha(c)) {
            if (this.xmlMode || this.sectionStart + 1 === this._index) ;
            else if (this.baseState !== 1 /* Text */) {
                if (c !== "=") {
                    // Parse as legacy entity, without allowing additional characters.
                    this.parseFixedEntity(legacy_json_1$1.default);
                }
            }
            else {
                this.parseLegacyEntity();
            }
            this._state = this.baseState;
            this._index--;
        }
    };
    Tokenizer.prototype.decodeNumericEntity = function (offset, base, strict) {
        var sectionStart = this.sectionStart + offset;
        if (sectionStart !== this._index) {
            // Parse entity
            var entity = this.buffer.substring(sectionStart, this._index);
            var parsed = parseInt(entity, base);
            this.emitPartial(decode_codepoint_1$1.default(parsed));
            this.sectionStart = strict ? this._index + 1 : this._index;
        }
        this._state = this.baseState;
    };
    Tokenizer.prototype.stateInNumericEntity = function (c) {
        if (c === ";") {
            this.decodeNumericEntity(2, 10, true);
        }
        else if (c < "0" || c > "9") {
            if (!this.xmlMode) {
                this.decodeNumericEntity(2, 10, false);
            }
            else {
                this._state = this.baseState;
            }
            this._index--;
        }
    };
    Tokenizer.prototype.stateInHexEntity = function (c) {
        if (c === ";") {
            this.decodeNumericEntity(3, 16, true);
        }
        else if ((c < "a" || c > "f") &&
            (c < "A" || c > "F") &&
            (c < "0" || c > "9")) {
            if (!this.xmlMode) {
                this.decodeNumericEntity(3, 16, false);
            }
            else {
                this._state = this.baseState;
            }
            this._index--;
        }
    };
    Tokenizer.prototype.cleanup = function () {
        if (this.sectionStart < 0) {
            this.buffer = "";
            this.bufferOffset += this._index;
            this._index = 0;
        }
        else if (this.running) {
            if (this._state === 1 /* Text */) {
                if (this.sectionStart !== this._index) {
                    this.cbs.ontext(this.buffer.substr(this.sectionStart));
                }
                this.buffer = "";
                this.bufferOffset += this._index;
                this._index = 0;
            }
            else if (this.sectionStart === this._index) {
                // The section just started
                this.buffer = "";
                this.bufferOffset += this._index;
                this._index = 0;
            }
            else {
                // Remove everything unnecessary
                this.buffer = this.buffer.substr(this.sectionStart);
                this._index -= this.sectionStart;
                this.bufferOffset += this.sectionStart;
            }
            this.sectionStart = 0;
        }
    };
    /**
     * Iterates through the buffer, calling the function corresponding to the current state.
     *
     * States that are more likely to be hit are higher up, as a performance improvement.
     */
    Tokenizer.prototype.parse = function () {
        while (this._index < this.buffer.length && this.running) {
            var c = this.buffer.charAt(this._index);
            if (this._state === 1 /* Text */) {
                this.stateText(c);
            }
            else if (this._state === 12 /* InAttributeValueDq */) {
                this.stateInAttributeValueDoubleQuotes(c);
            }
            else if (this._state === 9 /* InAttributeName */) {
                this.stateInAttributeName(c);
            }
            else if (this._state === 19 /* InComment */) {
                this.stateInComment(c);
            }
            else if (this._state === 20 /* InSpecialComment */) {
                this.stateInSpecialComment(c);
            }
            else if (this._state === 8 /* BeforeAttributeName */) {
                this.stateBeforeAttributeName(c);
            }
            else if (this._state === 3 /* InTagName */) {
                this.stateInTagName(c);
            }
            else if (this._state === 6 /* InClosingTagName */) {
                this.stateInClosingTagName(c);
            }
            else if (this._state === 2 /* BeforeTagName */) {
                this.stateBeforeTagName(c);
            }
            else if (this._state === 10 /* AfterAttributeName */) {
                this.stateAfterAttributeName(c);
            }
            else if (this._state === 13 /* InAttributeValueSq */) {
                this.stateInAttributeValueSingleQuotes(c);
            }
            else if (this._state === 11 /* BeforeAttributeValue */) {
                this.stateBeforeAttributeValue(c);
            }
            else if (this._state === 5 /* BeforeClosingTagName */) {
                this.stateBeforeClosingTagName(c);
            }
            else if (this._state === 7 /* AfterClosingTagName */) {
                this.stateAfterClosingTagName(c);
            }
            else if (this._state === 32 /* BeforeSpecialS */) {
                this.stateBeforeSpecialS(c);
            }
            else if (this._state === 21 /* AfterComment1 */) {
                this.stateAfterComment1(c);
            }
            else if (this._state === 14 /* InAttributeValueNq */) {
                this.stateInAttributeValueNoQuotes(c);
            }
            else if (this._state === 4 /* InSelfClosingTag */) {
                this.stateInSelfClosingTag(c);
            }
            else if (this._state === 16 /* InDeclaration */) {
                this.stateInDeclaration(c);
            }
            else if (this._state === 15 /* BeforeDeclaration */) {
                this.stateBeforeDeclaration(c);
            }
            else if (this._state === 22 /* AfterComment2 */) {
                this.stateAfterComment2(c);
            }
            else if (this._state === 18 /* BeforeComment */) {
                this.stateBeforeComment(c);
            }
            else if (this._state === 33 /* BeforeSpecialSEnd */) {
                this.stateBeforeSpecialSEnd(c);
            }
            else if (this._state === 53 /* BeforeSpecialTEnd */) {
                stateAfterSpecialTEnd(this, c);
            }
            else if (this._state === 39 /* AfterScript1 */) {
                stateAfterScript1(this, c);
            }
            else if (this._state === 40 /* AfterScript2 */) {
                stateAfterScript2(this, c);
            }
            else if (this._state === 41 /* AfterScript3 */) {
                stateAfterScript3(this, c);
            }
            else if (this._state === 34 /* BeforeScript1 */) {
                stateBeforeScript1(this, c);
            }
            else if (this._state === 35 /* BeforeScript2 */) {
                stateBeforeScript2(this, c);
            }
            else if (this._state === 36 /* BeforeScript3 */) {
                stateBeforeScript3(this, c);
            }
            else if (this._state === 37 /* BeforeScript4 */) {
                stateBeforeScript4(this, c);
            }
            else if (this._state === 38 /* BeforeScript5 */) {
                this.stateBeforeSpecialLast(c, 2 /* Script */);
            }
            else if (this._state === 42 /* AfterScript4 */) {
                stateAfterScript4(this, c);
            }
            else if (this._state === 43 /* AfterScript5 */) {
                this.stateAfterSpecialLast(c, 6);
            }
            else if (this._state === 44 /* BeforeStyle1 */) {
                stateBeforeStyle1(this, c);
            }
            else if (this._state === 29 /* InCdata */) {
                this.stateInCdata(c);
            }
            else if (this._state === 45 /* BeforeStyle2 */) {
                stateBeforeStyle2(this, c);
            }
            else if (this._state === 46 /* BeforeStyle3 */) {
                stateBeforeStyle3(this, c);
            }
            else if (this._state === 47 /* BeforeStyle4 */) {
                this.stateBeforeSpecialLast(c, 3 /* Style */);
            }
            else if (this._state === 48 /* AfterStyle1 */) {
                stateAfterStyle1(this, c);
            }
            else if (this._state === 49 /* AfterStyle2 */) {
                stateAfterStyle2(this, c);
            }
            else if (this._state === 50 /* AfterStyle3 */) {
                stateAfterStyle3(this, c);
            }
            else if (this._state === 51 /* AfterStyle4 */) {
                this.stateAfterSpecialLast(c, 5);
            }
            else if (this._state === 52 /* BeforeSpecialT */) {
                stateBeforeSpecialT(this, c);
            }
            else if (this._state === 54 /* BeforeTitle1 */) {
                stateBeforeTitle1(this, c);
            }
            else if (this._state === 55 /* BeforeTitle2 */) {
                stateBeforeTitle2(this, c);
            }
            else if (this._state === 56 /* BeforeTitle3 */) {
                stateBeforeTitle3(this, c);
            }
            else if (this._state === 57 /* BeforeTitle4 */) {
                this.stateBeforeSpecialLast(c, 4 /* Title */);
            }
            else if (this._state === 58 /* AfterTitle1 */) {
                stateAfterTitle1(this, c);
            }
            else if (this._state === 59 /* AfterTitle2 */) {
                stateAfterTitle2(this, c);
            }
            else if (this._state === 60 /* AfterTitle3 */) {
                stateAfterTitle3(this, c);
            }
            else if (this._state === 61 /* AfterTitle4 */) {
                this.stateAfterSpecialLast(c, 5);
            }
            else if (this._state === 17 /* InProcessingInstruction */) {
                this.stateInProcessingInstruction(c);
            }
            else if (this._state === 64 /* InNamedEntity */) {
                this.stateInNamedEntity(c);
            }
            else if (this._state === 23 /* BeforeCdata1 */) {
                stateBeforeCdata1(this, c);
            }
            else if (this._state === 62 /* BeforeEntity */) {
                stateBeforeEntity(this, c);
            }
            else if (this._state === 24 /* BeforeCdata2 */) {
                stateBeforeCdata2(this, c);
            }
            else if (this._state === 25 /* BeforeCdata3 */) {
                stateBeforeCdata3(this, c);
            }
            else if (this._state === 30 /* AfterCdata1 */) {
                this.stateAfterCdata1(c);
            }
            else if (this._state === 31 /* AfterCdata2 */) {
                this.stateAfterCdata2(c);
            }
            else if (this._state === 26 /* BeforeCdata4 */) {
                stateBeforeCdata4(this, c);
            }
            else if (this._state === 27 /* BeforeCdata5 */) {
                stateBeforeCdata5(this, c);
            }
            else if (this._state === 28 /* BeforeCdata6 */) {
                this.stateBeforeCdata6(c);
            }
            else if (this._state === 66 /* InHexEntity */) {
                this.stateInHexEntity(c);
            }
            else if (this._state === 65 /* InNumericEntity */) {
                this.stateInNumericEntity(c);
                // eslint-disable-next-line @typescript-eslint/no-unnecessary-condition
            }
            else if (this._state === 63 /* BeforeNumericEntity */) {
                stateBeforeNumericEntity(this, c);
            }
            else {
                this.cbs.onerror(Error("unknown _state"), this._state);
            }
            this._index++;
        }
        this.cleanup();
    };
    Tokenizer.prototype.finish = function () {
        // If there is remaining data, emit it in a reasonable way
        if (this.sectionStart < this._index) {
            this.handleTrailingData();
        }
        this.cbs.onend();
    };
    Tokenizer.prototype.handleTrailingData = function () {
        var data = this.buffer.substr(this.sectionStart);
        if (this._state === 29 /* InCdata */ ||
            this._state === 30 /* AfterCdata1 */ ||
            this._state === 31 /* AfterCdata2 */) {
            this.cbs.oncdata(data);
        }
        else if (this._state === 19 /* InComment */ ||
            this._state === 21 /* AfterComment1 */ ||
            this._state === 22 /* AfterComment2 */) {
            this.cbs.oncomment(data);
        }
        else if (this._state === 64 /* InNamedEntity */ && !this.xmlMode) {
            this.parseLegacyEntity();
            if (this.sectionStart < this._index) {
                this._state = this.baseState;
                this.handleTrailingData();
            }
        }
        else if (this._state === 65 /* InNumericEntity */ && !this.xmlMode) {
            this.decodeNumericEntity(2, 10, false);
            if (this.sectionStart < this._index) {
                this._state = this.baseState;
                this.handleTrailingData();
            }
        }
        else if (this._state === 66 /* InHexEntity */ && !this.xmlMode) {
            this.decodeNumericEntity(3, 16, false);
            if (this.sectionStart < this._index) {
                this._state = this.baseState;
                this.handleTrailingData();
            }
        }
        else if (this._state !== 3 /* InTagName */ &&
            this._state !== 8 /* BeforeAttributeName */ &&
            this._state !== 11 /* BeforeAttributeValue */ &&
            this._state !== 10 /* AfterAttributeName */ &&
            this._state !== 9 /* InAttributeName */ &&
            this._state !== 13 /* InAttributeValueSq */ &&
            this._state !== 12 /* InAttributeValueDq */ &&
            this._state !== 14 /* InAttributeValueNq */ &&
            this._state !== 6 /* InClosingTagName */) {
            this.cbs.ontext(data);
        }
        /*
         * Else, ignore remaining data
         * TODO add a way to remove current tag
         */
    };
    Tokenizer.prototype.getSection = function () {
        return this.buffer.substring(this.sectionStart, this._index);
    };
    Tokenizer.prototype.emitToken = function (name) {
        this.cbs[name](this.getSection());
        this.sectionStart = -1;
    };
    Tokenizer.prototype.emitPartial = function (value) {
        if (this.baseState !== 1 /* Text */) {
            this.cbs.onattribdata(value); // TODO implement the new event
        }
        else {
            this.cbs.ontext(value);
        }
    };
    return Tokenizer;
}());
Tokenizer$1.default = Tokenizer;

var __importDefault$5 = (commonjsGlobal && commonjsGlobal.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(Parser$1, "__esModule", { value: true });
Parser$1.Parser = void 0;
var Tokenizer_1 = __importDefault$5(Tokenizer$1);
var formTags = new Set([
    "input",
    "option",
    "optgroup",
    "select",
    "button",
    "datalist",
    "textarea",
]);
var pTag = new Set(["p"]);
var openImpliesClose = {
    tr: new Set(["tr", "th", "td"]),
    th: new Set(["th"]),
    td: new Set(["thead", "th", "td"]),
    body: new Set(["head", "link", "script"]),
    li: new Set(["li"]),
    p: pTag,
    h1: pTag,
    h2: pTag,
    h3: pTag,
    h4: pTag,
    h5: pTag,
    h6: pTag,
    select: formTags,
    input: formTags,
    output: formTags,
    button: formTags,
    datalist: formTags,
    textarea: formTags,
    option: new Set(["option"]),
    optgroup: new Set(["optgroup", "option"]),
    dd: new Set(["dt", "dd"]),
    dt: new Set(["dt", "dd"]),
    address: pTag,
    article: pTag,
    aside: pTag,
    blockquote: pTag,
    details: pTag,
    div: pTag,
    dl: pTag,
    fieldset: pTag,
    figcaption: pTag,
    figure: pTag,
    footer: pTag,
    form: pTag,
    header: pTag,
    hr: pTag,
    main: pTag,
    nav: pTag,
    ol: pTag,
    pre: pTag,
    section: pTag,
    table: pTag,
    ul: pTag,
    rt: new Set(["rt", "rp"]),
    rp: new Set(["rt", "rp"]),
    tbody: new Set(["thead", "tbody"]),
    tfoot: new Set(["thead", "tbody"]),
};
var voidElements = new Set([
    "area",
    "base",
    "basefont",
    "br",
    "col",
    "command",
    "embed",
    "frame",
    "hr",
    "img",
    "input",
    "isindex",
    "keygen",
    "link",
    "meta",
    "param",
    "source",
    "track",
    "wbr",
]);
var foreignContextElements = new Set(["math", "svg"]);
var htmlIntegrationElements = new Set([
    "mi",
    "mo",
    "mn",
    "ms",
    "mtext",
    "annotation-xml",
    "foreignObject",
    "desc",
    "title",
]);
var reNameEnd = /\s|\//;
var Parser = /** @class */ (function () {
    function Parser(cbs, options) {
        if (options === void 0) { options = {}; }
        var _a, _b, _c, _d, _e;
        /** The start index of the last event. */
        this.startIndex = 0;
        /** The end index of the last event. */
        this.endIndex = null;
        this.tagname = "";
        this.attribname = "";
        this.attribvalue = "";
        this.attribs = null;
        this.stack = [];
        this.foreignContext = [];
        this.options = options;
        this.cbs = cbs !== null && cbs !== void 0 ? cbs : {};
        this.lowerCaseTagNames = (_a = options.lowerCaseTags) !== null && _a !== void 0 ? _a : !options.xmlMode;
        this.lowerCaseAttributeNames =
            (_b = options.lowerCaseAttributeNames) !== null && _b !== void 0 ? _b : !options.xmlMode;
        this.tokenizer = new ((_c = options.Tokenizer) !== null && _c !== void 0 ? _c : Tokenizer_1.default)(this.options, this);
        (_e = (_d = this.cbs).onparserinit) === null || _e === void 0 ? void 0 : _e.call(_d, this);
    }
    Parser.prototype.updatePosition = function (initialOffset) {
        if (this.endIndex === null) {
            if (this.tokenizer.sectionStart <= initialOffset) {
                this.startIndex = 0;
            }
            else {
                this.startIndex = this.tokenizer.sectionStart - initialOffset;
            }
        }
        else {
            this.startIndex = this.endIndex + 1;
        }
        this.endIndex = this.tokenizer.getAbsoluteIndex();
    };
    // Tokenizer event handlers
    Parser.prototype.ontext = function (data) {
        var _a, _b;
        this.updatePosition(1);
        this.endIndex--;
        (_b = (_a = this.cbs).ontext) === null || _b === void 0 ? void 0 : _b.call(_a, data);
    };
    Parser.prototype.onopentagname = function (name) {
        var _a, _b;
        if (this.lowerCaseTagNames) {
            name = name.toLowerCase();
        }
        this.tagname = name;
        if (!this.options.xmlMode &&
            Object.prototype.hasOwnProperty.call(openImpliesClose, name)) {
            var el = void 0;
            while (this.stack.length > 0 &&
                openImpliesClose[name].has((el = this.stack[this.stack.length - 1]))) {
                this.onclosetag(el);
            }
        }
        if (this.options.xmlMode || !voidElements.has(name)) {
            this.stack.push(name);
            if (foreignContextElements.has(name)) {
                this.foreignContext.push(true);
            }
            else if (htmlIntegrationElements.has(name)) {
                this.foreignContext.push(false);
            }
        }
        (_b = (_a = this.cbs).onopentagname) === null || _b === void 0 ? void 0 : _b.call(_a, name);
        if (this.cbs.onopentag)
            this.attribs = {};
    };
    Parser.prototype.onopentagend = function () {
        var _a, _b;
        this.updatePosition(1);
        if (this.attribs) {
            (_b = (_a = this.cbs).onopentag) === null || _b === void 0 ? void 0 : _b.call(_a, this.tagname, this.attribs);
            this.attribs = null;
        }
        if (!this.options.xmlMode &&
            this.cbs.onclosetag &&
            voidElements.has(this.tagname)) {
            this.cbs.onclosetag(this.tagname);
        }
        this.tagname = "";
    };
    Parser.prototype.onclosetag = function (name) {
        this.updatePosition(1);
        if (this.lowerCaseTagNames) {
            name = name.toLowerCase();
        }
        if (foreignContextElements.has(name) ||
            htmlIntegrationElements.has(name)) {
            this.foreignContext.pop();
        }
        if (this.stack.length &&
            (this.options.xmlMode || !voidElements.has(name))) {
            var pos = this.stack.lastIndexOf(name);
            if (pos !== -1) {
                if (this.cbs.onclosetag) {
                    pos = this.stack.length - pos;
                    while (pos--) {
                        // We know the stack has sufficient elements.
                        this.cbs.onclosetag(this.stack.pop());
                    }
                }
                else
                    this.stack.length = pos;
            }
            else if (name === "p" && !this.options.xmlMode) {
                this.onopentagname(name);
                this.closeCurrentTag();
            }
        }
        else if (!this.options.xmlMode && (name === "br" || name === "p")) {
            this.onopentagname(name);
            this.closeCurrentTag();
        }
    };
    Parser.prototype.onselfclosingtag = function () {
        if (this.options.xmlMode ||
            this.options.recognizeSelfClosing ||
            this.foreignContext[this.foreignContext.length - 1]) {
            this.closeCurrentTag();
        }
        else {
            this.onopentagend();
        }
    };
    Parser.prototype.closeCurrentTag = function () {
        var _a, _b;
        var name = this.tagname;
        this.onopentagend();
        /*
         * Self-closing tags will be on the top of the stack
         * (cheaper check than in onclosetag)
         */
        if (this.stack[this.stack.length - 1] === name) {
            (_b = (_a = this.cbs).onclosetag) === null || _b === void 0 ? void 0 : _b.call(_a, name);
            this.stack.pop();
        }
    };
    Parser.prototype.onattribname = function (name) {
        if (this.lowerCaseAttributeNames) {
            name = name.toLowerCase();
        }
        this.attribname = name;
    };
    Parser.prototype.onattribdata = function (value) {
        this.attribvalue += value;
    };
    Parser.prototype.onattribend = function (quote) {
        var _a, _b;
        (_b = (_a = this.cbs).onattribute) === null || _b === void 0 ? void 0 : _b.call(_a, this.attribname, this.attribvalue, quote);
        if (this.attribs &&
            !Object.prototype.hasOwnProperty.call(this.attribs, this.attribname)) {
            this.attribs[this.attribname] = this.attribvalue;
        }
        this.attribname = "";
        this.attribvalue = "";
    };
    Parser.prototype.getInstructionName = function (value) {
        var idx = value.search(reNameEnd);
        var name = idx < 0 ? value : value.substr(0, idx);
        if (this.lowerCaseTagNames) {
            name = name.toLowerCase();
        }
        return name;
    };
    Parser.prototype.ondeclaration = function (value) {
        if (this.cbs.onprocessinginstruction) {
            var name_1 = this.getInstructionName(value);
            this.cbs.onprocessinginstruction("!" + name_1, "!" + value);
        }
    };
    Parser.prototype.onprocessinginstruction = function (value) {
        if (this.cbs.onprocessinginstruction) {
            var name_2 = this.getInstructionName(value);
            this.cbs.onprocessinginstruction("?" + name_2, "?" + value);
        }
    };
    Parser.prototype.oncomment = function (value) {
        var _a, _b, _c, _d;
        this.updatePosition(4);
        (_b = (_a = this.cbs).oncomment) === null || _b === void 0 ? void 0 : _b.call(_a, value);
        (_d = (_c = this.cbs).oncommentend) === null || _d === void 0 ? void 0 : _d.call(_c);
    };
    Parser.prototype.oncdata = function (value) {
        var _a, _b, _c, _d, _e, _f;
        this.updatePosition(1);
        if (this.options.xmlMode || this.options.recognizeCDATA) {
            (_b = (_a = this.cbs).oncdatastart) === null || _b === void 0 ? void 0 : _b.call(_a);
            (_d = (_c = this.cbs).ontext) === null || _d === void 0 ? void 0 : _d.call(_c, value);
            (_f = (_e = this.cbs).oncdataend) === null || _f === void 0 ? void 0 : _f.call(_e);
        }
        else {
            this.oncomment("[CDATA[" + value + "]]");
        }
    };
    Parser.prototype.onerror = function (err) {
        var _a, _b;
        (_b = (_a = this.cbs).onerror) === null || _b === void 0 ? void 0 : _b.call(_a, err);
    };
    Parser.prototype.onend = function () {
        var _a, _b;
        if (this.cbs.onclosetag) {
            for (var i = this.stack.length; i > 0; this.cbs.onclosetag(this.stack[--i]))
                ;
        }
        (_b = (_a = this.cbs).onend) === null || _b === void 0 ? void 0 : _b.call(_a);
    };
    /**
     * Resets the parser to a blank state, ready to parse a new HTML document
     */
    Parser.prototype.reset = function () {
        var _a, _b, _c, _d;
        (_b = (_a = this.cbs).onreset) === null || _b === void 0 ? void 0 : _b.call(_a);
        this.tokenizer.reset();
        this.tagname = "";
        this.attribname = "";
        this.attribs = null;
        this.stack = [];
        (_d = (_c = this.cbs).onparserinit) === null || _d === void 0 ? void 0 : _d.call(_c, this);
    };
    /**
     * Resets the parser, then parses a complete document and
     * pushes it to the handler.
     *
     * @param data Document to parse.
     */
    Parser.prototype.parseComplete = function (data) {
        this.reset();
        this.end(data);
    };
    /**
     * Parses a chunk of data and calls the corresponding callbacks.
     *
     * @param chunk Chunk to parse.
     */
    Parser.prototype.write = function (chunk) {
        this.tokenizer.write(chunk);
    };
    /**
     * Parses the end of the buffer and clears the stack, calls onend.
     *
     * @param chunk Optional final chunk to parse.
     */
    Parser.prototype.end = function (chunk) {
        this.tokenizer.end(chunk);
    };
    /**
     * Pauses parsing. The parser won't emit events until `resume` is called.
     */
    Parser.prototype.pause = function () {
        this.tokenizer.pause();
    };
    /**
     * Resumes parsing after `pause` was called.
     */
    Parser.prototype.resume = function () {
        this.tokenizer.resume();
    };
    /**
     * Alias of `write`, for backwards compatibility.
     *
     * @param chunk Chunk to parse.
     * @deprecated
     */
    Parser.prototype.parseChunk = function (chunk) {
        this.write(chunk);
    };
    /**
     * Alias of `end`, for backwards compatibility.
     *
     * @param chunk Optional final chunk to parse.
     * @deprecated
     */
    Parser.prototype.done = function (chunk) {
        this.end(chunk);
    };
    return Parser;
}());
Parser$1.Parser = Parser;

var FeedHandler$1 = {};

var lib$2 = {};

var stringify = {};

var lib$1 = {};

var lib = {};

var decode = {};

var Aacute$1 = "Á";
var aacute$1 = "á";
var Abreve = "Ă";
var abreve = "ă";
var ac = "∾";
var acd = "∿";
var acE = "∾̳";
var Acirc$1 = "Â";
var acirc$1 = "â";
var acute$1 = "´";
var Acy = "А";
var acy = "а";
var AElig$1 = "Æ";
var aelig$1 = "æ";
var af = "⁡";
var Afr = "𝔄";
var afr = "𝔞";
var Agrave$1 = "À";
var agrave$1 = "à";
var alefsym = "ℵ";
var aleph = "ℵ";
var Alpha = "Α";
var alpha = "α";
var Amacr = "Ā";
var amacr = "ā";
var amalg = "⨿";
var amp$2 = "&";
var AMP$1 = "&";
var andand = "⩕";
var And = "⩓";
var and = "∧";
var andd = "⩜";
var andslope = "⩘";
var andv = "⩚";
var ang = "∠";
var ange = "⦤";
var angle = "∠";
var angmsdaa = "⦨";
var angmsdab = "⦩";
var angmsdac = "⦪";
var angmsdad = "⦫";
var angmsdae = "⦬";
var angmsdaf = "⦭";
var angmsdag = "⦮";
var angmsdah = "⦯";
var angmsd = "∡";
var angrt = "∟";
var angrtvb = "⊾";
var angrtvbd = "⦝";
var angsph = "∢";
var angst = "Å";
var angzarr = "⍼";
var Aogon = "Ą";
var aogon = "ą";
var Aopf = "𝔸";
var aopf = "𝕒";
var apacir = "⩯";
var ap = "≈";
var apE = "⩰";
var ape = "≊";
var apid = "≋";
var apos$1 = "'";
var ApplyFunction = "⁡";
var approx = "≈";
var approxeq = "≊";
var Aring$1 = "Å";
var aring$1 = "å";
var Ascr = "𝒜";
var ascr = "𝒶";
var Assign = "≔";
var ast = "*";
var asymp = "≈";
var asympeq = "≍";
var Atilde$1 = "Ã";
var atilde$1 = "ã";
var Auml$1 = "Ä";
var auml$1 = "ä";
var awconint = "∳";
var awint = "⨑";
var backcong = "≌";
var backepsilon = "϶";
var backprime = "‵";
var backsim = "∽";
var backsimeq = "⋍";
var Backslash = "∖";
var Barv = "⫧";
var barvee = "⊽";
var barwed = "⌅";
var Barwed = "⌆";
var barwedge = "⌅";
var bbrk = "⎵";
var bbrktbrk = "⎶";
var bcong = "≌";
var Bcy = "Б";
var bcy = "б";
var bdquo = "„";
var becaus = "∵";
var because = "∵";
var Because = "∵";
var bemptyv = "⦰";
var bepsi = "϶";
var bernou = "ℬ";
var Bernoullis = "ℬ";
var Beta = "Β";
var beta = "β";
var beth = "ℶ";
var between = "≬";
var Bfr = "𝔅";
var bfr = "𝔟";
var bigcap = "⋂";
var bigcirc = "◯";
var bigcup = "⋃";
var bigodot = "⨀";
var bigoplus = "⨁";
var bigotimes = "⨂";
var bigsqcup = "⨆";
var bigstar = "★";
var bigtriangledown = "▽";
var bigtriangleup = "△";
var biguplus = "⨄";
var bigvee = "⋁";
var bigwedge = "⋀";
var bkarow = "⤍";
var blacklozenge = "⧫";
var blacksquare = "▪";
var blacktriangle = "▴";
var blacktriangledown = "▾";
var blacktriangleleft = "◂";
var blacktriangleright = "▸";
var blank = "␣";
var blk12 = "▒";
var blk14 = "░";
var blk34 = "▓";
var block = "█";
var bne = "=⃥";
var bnequiv = "≡⃥";
var bNot = "⫭";
var bnot = "⌐";
var Bopf = "𝔹";
var bopf = "𝕓";
var bot = "⊥";
var bottom = "⊥";
var bowtie = "⋈";
var boxbox = "⧉";
var boxdl = "┐";
var boxdL = "╕";
var boxDl = "╖";
var boxDL = "╗";
var boxdr = "┌";
var boxdR = "╒";
var boxDr = "╓";
var boxDR = "╔";
var boxh = "─";
var boxH = "═";
var boxhd = "┬";
var boxHd = "╤";
var boxhD = "╥";
var boxHD = "╦";
var boxhu = "┴";
var boxHu = "╧";
var boxhU = "╨";
var boxHU = "╩";
var boxminus = "⊟";
var boxplus = "⊞";
var boxtimes = "⊠";
var boxul = "┘";
var boxuL = "╛";
var boxUl = "╜";
var boxUL = "╝";
var boxur = "└";
var boxuR = "╘";
var boxUr = "╙";
var boxUR = "╚";
var boxv = "│";
var boxV = "║";
var boxvh = "┼";
var boxvH = "╪";
var boxVh = "╫";
var boxVH = "╬";
var boxvl = "┤";
var boxvL = "╡";
var boxVl = "╢";
var boxVL = "╣";
var boxvr = "├";
var boxvR = "╞";
var boxVr = "╟";
var boxVR = "╠";
var bprime = "‵";
var breve = "˘";
var Breve = "˘";
var brvbar$1 = "¦";
var bscr = "𝒷";
var Bscr = "ℬ";
var bsemi = "⁏";
var bsim = "∽";
var bsime = "⋍";
var bsolb = "⧅";
var bsol = "\\";
var bsolhsub = "⟈";
var bull = "•";
var bullet = "•";
var bump = "≎";
var bumpE = "⪮";
var bumpe = "≏";
var Bumpeq = "≎";
var bumpeq = "≏";
var Cacute = "Ć";
var cacute = "ć";
var capand = "⩄";
var capbrcup = "⩉";
var capcap = "⩋";
var cap = "∩";
var Cap = "⋒";
var capcup = "⩇";
var capdot = "⩀";
var CapitalDifferentialD = "ⅅ";
var caps = "∩︀";
var caret$1 = "⁁";
var caron = "ˇ";
var Cayleys = "ℭ";
var ccaps = "⩍";
var Ccaron = "Č";
var ccaron = "č";
var Ccedil$1 = "Ç";
var ccedil$1 = "ç";
var Ccirc = "Ĉ";
var ccirc = "ĉ";
var Cconint = "∰";
var ccups = "⩌";
var ccupssm = "⩐";
var Cdot = "Ċ";
var cdot = "ċ";
var cedil$1 = "¸";
var Cedilla = "¸";
var cemptyv = "⦲";
var cent$1 = "¢";
var centerdot = "·";
var CenterDot = "·";
var cfr = "𝔠";
var Cfr = "ℭ";
var CHcy = "Ч";
var chcy = "ч";
var check = "✓";
var checkmark = "✓";
var Chi = "Χ";
var chi = "χ";
var circ = "ˆ";
var circeq = "≗";
var circlearrowleft = "↺";
var circlearrowright = "↻";
var circledast = "⊛";
var circledcirc = "⊚";
var circleddash = "⊝";
var CircleDot = "⊙";
var circledR = "®";
var circledS = "Ⓢ";
var CircleMinus = "⊖";
var CirclePlus = "⊕";
var CircleTimes = "⊗";
var cir = "○";
var cirE = "⧃";
var cire = "≗";
var cirfnint = "⨐";
var cirmid = "⫯";
var cirscir = "⧂";
var ClockwiseContourIntegral = "∲";
var CloseCurlyDoubleQuote = "”";
var CloseCurlyQuote = "’";
var clubs = "♣";
var clubsuit = "♣";
var colon = ":";
var Colon = "∷";
var Colone = "⩴";
var colone = "≔";
var coloneq = "≔";
var comma$1 = ",";
var commat = "@";
var comp = "∁";
var compfn = "∘";
var complement = "∁";
var complexes = "ℂ";
var cong = "≅";
var congdot = "⩭";
var Congruent = "≡";
var conint = "∮";
var Conint = "∯";
var ContourIntegral = "∮";
var copf = "𝕔";
var Copf = "ℂ";
var coprod = "∐";
var Coproduct = "∐";
var copy$1 = "©";
var COPY$1 = "©";
var copysr = "℗";
var CounterClockwiseContourIntegral = "∳";
var crarr = "↵";
var cross = "✗";
var Cross = "⨯";
var Cscr = "𝒞";
var cscr = "𝒸";
var csub = "⫏";
var csube = "⫑";
var csup = "⫐";
var csupe = "⫒";
var ctdot = "⋯";
var cudarrl = "⤸";
var cudarrr = "⤵";
var cuepr = "⋞";
var cuesc = "⋟";
var cularr = "↶";
var cularrp = "⤽";
var cupbrcap = "⩈";
var cupcap = "⩆";
var CupCap = "≍";
var cup = "∪";
var Cup = "⋓";
var cupcup = "⩊";
var cupdot = "⊍";
var cupor = "⩅";
var cups = "∪︀";
var curarr = "↷";
var curarrm = "⤼";
var curlyeqprec = "⋞";
var curlyeqsucc = "⋟";
var curlyvee = "⋎";
var curlywedge = "⋏";
var curren$1 = "¤";
var curvearrowleft = "↶";
var curvearrowright = "↷";
var cuvee = "⋎";
var cuwed = "⋏";
var cwconint = "∲";
var cwint = "∱";
var cylcty = "⌭";
var dagger = "†";
var Dagger = "‡";
var daleth = "ℸ";
var darr = "↓";
var Darr = "↡";
var dArr = "⇓";
var dash = "‐";
var Dashv = "⫤";
var dashv = "⊣";
var dbkarow = "⤏";
var dblac = "˝";
var Dcaron = "Ď";
var dcaron = "ď";
var Dcy = "Д";
var dcy = "д";
var ddagger = "‡";
var ddarr = "⇊";
var DD = "ⅅ";
var dd = "ⅆ";
var DDotrahd = "⤑";
var ddotseq = "⩷";
var deg$1 = "°";
var Del = "∇";
var Delta = "Δ";
var delta = "δ";
var demptyv = "⦱";
var dfisht = "⥿";
var Dfr = "𝔇";
var dfr = "𝔡";
var dHar = "⥥";
var dharl = "⇃";
var dharr = "⇂";
var DiacriticalAcute = "´";
var DiacriticalDot = "˙";
var DiacriticalDoubleAcute = "˝";
var DiacriticalGrave = "`";
var DiacriticalTilde = "˜";
var diam = "⋄";
var diamond = "⋄";
var Diamond = "⋄";
var diamondsuit = "♦";
var diams = "♦";
var die = "¨";
var DifferentialD = "ⅆ";
var digamma = "ϝ";
var disin = "⋲";
var div = "÷";
var divide$1 = "÷";
var divideontimes = "⋇";
var divonx = "⋇";
var DJcy = "Ђ";
var djcy = "ђ";
var dlcorn = "⌞";
var dlcrop = "⌍";
var dollar$1 = "$";
var Dopf = "𝔻";
var dopf = "𝕕";
var Dot = "¨";
var dot = "˙";
var DotDot = "⃜";
var doteq = "≐";
var doteqdot = "≑";
var DotEqual = "≐";
var dotminus = "∸";
var dotplus = "∔";
var dotsquare = "⊡";
var doublebarwedge = "⌆";
var DoubleContourIntegral = "∯";
var DoubleDot = "¨";
var DoubleDownArrow = "⇓";
var DoubleLeftArrow = "⇐";
var DoubleLeftRightArrow = "⇔";
var DoubleLeftTee = "⫤";
var DoubleLongLeftArrow = "⟸";
var DoubleLongLeftRightArrow = "⟺";
var DoubleLongRightArrow = "⟹";
var DoubleRightArrow = "⇒";
var DoubleRightTee = "⊨";
var DoubleUpArrow = "⇑";
var DoubleUpDownArrow = "⇕";
var DoubleVerticalBar = "∥";
var DownArrowBar = "⤓";
var downarrow = "↓";
var DownArrow = "↓";
var Downarrow = "⇓";
var DownArrowUpArrow = "⇵";
var DownBreve = "̑";
var downdownarrows = "⇊";
var downharpoonleft = "⇃";
var downharpoonright = "⇂";
var DownLeftRightVector = "⥐";
var DownLeftTeeVector = "⥞";
var DownLeftVectorBar = "⥖";
var DownLeftVector = "↽";
var DownRightTeeVector = "⥟";
var DownRightVectorBar = "⥗";
var DownRightVector = "⇁";
var DownTeeArrow = "↧";
var DownTee = "⊤";
var drbkarow = "⤐";
var drcorn = "⌟";
var drcrop = "⌌";
var Dscr = "𝒟";
var dscr = "𝒹";
var DScy = "Ѕ";
var dscy = "ѕ";
var dsol = "⧶";
var Dstrok = "Đ";
var dstrok = "đ";
var dtdot = "⋱";
var dtri = "▿";
var dtrif = "▾";
var duarr = "⇵";
var duhar = "⥯";
var dwangle = "⦦";
var DZcy = "Џ";
var dzcy = "џ";
var dzigrarr = "⟿";
var Eacute$1 = "É";
var eacute$1 = "é";
var easter = "⩮";
var Ecaron = "Ě";
var ecaron = "ě";
var Ecirc$1 = "Ê";
var ecirc$1 = "ê";
var ecir = "≖";
var ecolon = "≕";
var Ecy = "Э";
var ecy = "э";
var eDDot = "⩷";
var Edot = "Ė";
var edot = "ė";
var eDot = "≑";
var ee = "ⅇ";
var efDot = "≒";
var Efr = "𝔈";
var efr = "𝔢";
var eg = "⪚";
var Egrave$1 = "È";
var egrave$1 = "è";
var egs = "⪖";
var egsdot = "⪘";
var el = "⪙";
var Element = "∈";
var elinters = "⏧";
var ell = "ℓ";
var els = "⪕";
var elsdot = "⪗";
var Emacr = "Ē";
var emacr = "ē";
var empty = "∅";
var emptyset = "∅";
var EmptySmallSquare = "◻";
var emptyv = "∅";
var EmptyVerySmallSquare = "▫";
var emsp13 = " ";
var emsp14 = " ";
var emsp = " ";
var ENG = "Ŋ";
var eng = "ŋ";
var ensp = " ";
var Eogon = "Ę";
var eogon = "ę";
var Eopf = "𝔼";
var eopf = "𝕖";
var epar = "⋕";
var eparsl = "⧣";
var eplus = "⩱";
var epsi = "ε";
var Epsilon = "Ε";
var epsilon = "ε";
var epsiv = "ϵ";
var eqcirc = "≖";
var eqcolon = "≕";
var eqsim = "≂";
var eqslantgtr = "⪖";
var eqslantless = "⪕";
var Equal = "⩵";
var equals = "=";
var EqualTilde = "≂";
var equest = "≟";
var Equilibrium = "⇌";
var equiv = "≡";
var equivDD = "⩸";
var eqvparsl = "⧥";
var erarr = "⥱";
var erDot = "≓";
var escr = "ℯ";
var Escr = "ℰ";
var esdot = "≐";
var Esim = "⩳";
var esim = "≂";
var Eta = "Η";
var eta = "η";
var ETH$1 = "Ð";
var eth$1 = "ð";
var Euml$1 = "Ë";
var euml$1 = "ë";
var euro = "€";
var excl = "!";
var exist = "∃";
var Exists = "∃";
var expectation = "ℰ";
var exponentiale = "ⅇ";
var ExponentialE = "ⅇ";
var fallingdotseq = "≒";
var Fcy = "Ф";
var fcy = "ф";
var female = "♀";
var ffilig = "ﬃ";
var fflig = "ﬀ";
var ffllig = "ﬄ";
var Ffr = "𝔉";
var ffr = "𝔣";
var filig = "ﬁ";
var FilledSmallSquare = "◼";
var FilledVerySmallSquare = "▪";
var fjlig = "fj";
var flat = "♭";
var fllig = "ﬂ";
var fltns = "▱";
var fnof = "ƒ";
var Fopf = "𝔽";
var fopf = "𝕗";
var forall = "∀";
var ForAll = "∀";
var fork = "⋔";
var forkv = "⫙";
var Fouriertrf = "ℱ";
var fpartint = "⨍";
var frac12$1 = "½";
var frac13 = "⅓";
var frac14$1 = "¼";
var frac15 = "⅕";
var frac16 = "⅙";
var frac18 = "⅛";
var frac23 = "⅔";
var frac25 = "⅖";
var frac34$1 = "¾";
var frac35 = "⅗";
var frac38 = "⅜";
var frac45 = "⅘";
var frac56 = "⅚";
var frac58 = "⅝";
var frac78 = "⅞";
var frasl = "⁄";
var frown = "⌢";
var fscr = "𝒻";
var Fscr = "ℱ";
var gacute = "ǵ";
var Gamma = "Γ";
var gamma = "γ";
var Gammad = "Ϝ";
var gammad = "ϝ";
var gap = "⪆";
var Gbreve = "Ğ";
var gbreve = "ğ";
var Gcedil = "Ģ";
var Gcirc = "Ĝ";
var gcirc = "ĝ";
var Gcy = "Г";
var gcy = "г";
var Gdot = "Ġ";
var gdot = "ġ";
var ge = "≥";
var gE = "≧";
var gEl = "⪌";
var gel = "⋛";
var geq = "≥";
var geqq = "≧";
var geqslant = "⩾";
var gescc = "⪩";
var ges = "⩾";
var gesdot = "⪀";
var gesdoto = "⪂";
var gesdotol = "⪄";
var gesl = "⋛︀";
var gesles = "⪔";
var Gfr = "𝔊";
var gfr = "𝔤";
var gg = "≫";
var Gg = "⋙";
var ggg = "⋙";
var gimel = "ℷ";
var GJcy = "Ѓ";
var gjcy = "ѓ";
var gla = "⪥";
var gl = "≷";
var glE = "⪒";
var glj = "⪤";
var gnap = "⪊";
var gnapprox = "⪊";
var gne = "⪈";
var gnE = "≩";
var gneq = "⪈";
var gneqq = "≩";
var gnsim = "⋧";
var Gopf = "𝔾";
var gopf = "𝕘";
var grave = "`";
var GreaterEqual = "≥";
var GreaterEqualLess = "⋛";
var GreaterFullEqual = "≧";
var GreaterGreater = "⪢";
var GreaterLess = "≷";
var GreaterSlantEqual = "⩾";
var GreaterTilde = "≳";
var Gscr = "𝒢";
var gscr = "ℊ";
var gsim = "≳";
var gsime = "⪎";
var gsiml = "⪐";
var gtcc = "⪧";
var gtcir = "⩺";
var gt$3 = ">";
var GT$1 = ">";
var Gt = "≫";
var gtdot = "⋗";
var gtlPar = "⦕";
var gtquest = "⩼";
var gtrapprox = "⪆";
var gtrarr = "⥸";
var gtrdot = "⋗";
var gtreqless = "⋛";
var gtreqqless = "⪌";
var gtrless = "≷";
var gtrsim = "≳";
var gvertneqq = "≩︀";
var gvnE = "≩︀";
var Hacek = "ˇ";
var hairsp = " ";
var half = "½";
var hamilt = "ℋ";
var HARDcy = "Ъ";
var hardcy = "ъ";
var harrcir = "⥈";
var harr = "↔";
var hArr = "⇔";
var harrw = "↭";
var Hat = "^";
var hbar = "ℏ";
var Hcirc = "Ĥ";
var hcirc = "ĥ";
var hearts = "♥";
var heartsuit = "♥";
var hellip = "…";
var hercon = "⊹";
var hfr = "𝔥";
var Hfr = "ℌ";
var HilbertSpace = "ℋ";
var hksearow = "⤥";
var hkswarow = "⤦";
var hoarr = "⇿";
var homtht = "∻";
var hookleftarrow = "↩";
var hookrightarrow = "↪";
var hopf = "𝕙";
var Hopf = "ℍ";
var horbar = "―";
var HorizontalLine = "─";
var hscr = "𝒽";
var Hscr = "ℋ";
var hslash = "ℏ";
var Hstrok = "Ħ";
var hstrok = "ħ";
var HumpDownHump = "≎";
var HumpEqual = "≏";
var hybull = "⁃";
var hyphen = "‐";
var Iacute$1 = "Í";
var iacute$1 = "í";
var ic = "⁣";
var Icirc$1 = "Î";
var icirc$1 = "î";
var Icy = "И";
var icy = "и";
var Idot = "İ";
var IEcy = "Е";
var iecy = "е";
var iexcl$1 = "¡";
var iff = "⇔";
var ifr = "𝔦";
var Ifr = "ℑ";
var Igrave$1 = "Ì";
var igrave$1 = "ì";
var ii = "ⅈ";
var iiiint = "⨌";
var iiint = "∭";
var iinfin = "⧜";
var iiota = "℩";
var IJlig = "Ĳ";
var ijlig = "ĳ";
var Imacr = "Ī";
var imacr = "ī";
var image = "ℑ";
var ImaginaryI = "ⅈ";
var imagline = "ℐ";
var imagpart = "ℑ";
var imath = "ı";
var Im = "ℑ";
var imof = "⊷";
var imped = "Ƶ";
var Implies = "⇒";
var incare = "℅";
var infin = "∞";
var infintie = "⧝";
var inodot = "ı";
var intcal = "⊺";
var int = "∫";
var Int = "∬";
var integers = "ℤ";
var Integral = "∫";
var intercal = "⊺";
var Intersection = "⋂";
var intlarhk = "⨗";
var intprod = "⨼";
var InvisibleComma = "⁣";
var InvisibleTimes = "⁢";
var IOcy = "Ё";
var iocy = "ё";
var Iogon = "Į";
var iogon = "į";
var Iopf = "𝕀";
var iopf = "𝕚";
var Iota = "Ι";
var iota = "ι";
var iprod = "⨼";
var iquest$1 = "¿";
var iscr = "𝒾";
var Iscr = "ℐ";
var isin = "∈";
var isindot = "⋵";
var isinE = "⋹";
var isins = "⋴";
var isinsv = "⋳";
var isinv = "∈";
var it = "⁢";
var Itilde = "Ĩ";
var itilde = "ĩ";
var Iukcy = "І";
var iukcy = "і";
var Iuml$1 = "Ï";
var iuml$1 = "ï";
var Jcirc = "Ĵ";
var jcirc = "ĵ";
var Jcy = "Й";
var jcy = "й";
var Jfr = "𝔍";
var jfr = "𝔧";
var jmath = "ȷ";
var Jopf = "𝕁";
var jopf = "𝕛";
var Jscr = "𝒥";
var jscr = "𝒿";
var Jsercy = "Ј";
var jsercy = "ј";
var Jukcy = "Є";
var jukcy = "є";
var Kappa = "Κ";
var kappa = "κ";
var kappav = "ϰ";
var Kcedil = "Ķ";
var kcedil = "ķ";
var Kcy = "К";
var kcy = "к";
var Kfr = "𝔎";
var kfr = "𝔨";
var kgreen = "ĸ";
var KHcy = "Х";
var khcy = "х";
var KJcy = "Ќ";
var kjcy = "ќ";
var Kopf = "𝕂";
var kopf = "𝕜";
var Kscr = "𝒦";
var kscr = "𝓀";
var lAarr = "⇚";
var Lacute = "Ĺ";
var lacute = "ĺ";
var laemptyv = "⦴";
var lagran = "ℒ";
var Lambda = "Λ";
var lambda = "λ";
var lang = "⟨";
var Lang = "⟪";
var langd = "⦑";
var langle = "⟨";
var lap = "⪅";
var Laplacetrf = "ℒ";
var laquo$1 = "«";
var larrb = "⇤";
var larrbfs = "⤟";
var larr = "←";
var Larr = "↞";
var lArr = "⇐";
var larrfs = "⤝";
var larrhk = "↩";
var larrlp = "↫";
var larrpl = "⤹";
var larrsim = "⥳";
var larrtl = "↢";
var latail = "⤙";
var lAtail = "⤛";
var lat = "⪫";
var late = "⪭";
var lates = "⪭︀";
var lbarr = "⤌";
var lBarr = "⤎";
var lbbrk = "❲";
var lbrace = "{";
var lbrack = "[";
var lbrke = "⦋";
var lbrksld = "⦏";
var lbrkslu = "⦍";
var Lcaron = "Ľ";
var lcaron = "ľ";
var Lcedil = "Ļ";
var lcedil = "ļ";
var lceil = "⌈";
var lcub = "{";
var Lcy = "Л";
var lcy = "л";
var ldca = "⤶";
var ldquo = "“";
var ldquor = "„";
var ldrdhar = "⥧";
var ldrushar = "⥋";
var ldsh = "↲";
var le = "≤";
var lE = "≦";
var LeftAngleBracket = "⟨";
var LeftArrowBar = "⇤";
var leftarrow = "←";
var LeftArrow = "←";
var Leftarrow = "⇐";
var LeftArrowRightArrow = "⇆";
var leftarrowtail = "↢";
var LeftCeiling = "⌈";
var LeftDoubleBracket = "⟦";
var LeftDownTeeVector = "⥡";
var LeftDownVectorBar = "⥙";
var LeftDownVector = "⇃";
var LeftFloor = "⌊";
var leftharpoondown = "↽";
var leftharpoonup = "↼";
var leftleftarrows = "⇇";
var leftrightarrow = "↔";
var LeftRightArrow = "↔";
var Leftrightarrow = "⇔";
var leftrightarrows = "⇆";
var leftrightharpoons = "⇋";
var leftrightsquigarrow = "↭";
var LeftRightVector = "⥎";
var LeftTeeArrow = "↤";
var LeftTee = "⊣";
var LeftTeeVector = "⥚";
var leftthreetimes = "⋋";
var LeftTriangleBar = "⧏";
var LeftTriangle = "⊲";
var LeftTriangleEqual = "⊴";
var LeftUpDownVector = "⥑";
var LeftUpTeeVector = "⥠";
var LeftUpVectorBar = "⥘";
var LeftUpVector = "↿";
var LeftVectorBar = "⥒";
var LeftVector = "↼";
var lEg = "⪋";
var leg = "⋚";
var leq = "≤";
var leqq = "≦";
var leqslant = "⩽";
var lescc = "⪨";
var les = "⩽";
var lesdot = "⩿";
var lesdoto = "⪁";
var lesdotor = "⪃";
var lesg = "⋚︀";
var lesges = "⪓";
var lessapprox = "⪅";
var lessdot = "⋖";
var lesseqgtr = "⋚";
var lesseqqgtr = "⪋";
var LessEqualGreater = "⋚";
var LessFullEqual = "≦";
var LessGreater = "≶";
var lessgtr = "≶";
var LessLess = "⪡";
var lesssim = "≲";
var LessSlantEqual = "⩽";
var LessTilde = "≲";
var lfisht = "⥼";
var lfloor = "⌊";
var Lfr = "𝔏";
var lfr = "𝔩";
var lg = "≶";
var lgE = "⪑";
var lHar = "⥢";
var lhard = "↽";
var lharu = "↼";
var lharul = "⥪";
var lhblk = "▄";
var LJcy = "Љ";
var ljcy = "љ";
var llarr = "⇇";
var ll = "≪";
var Ll = "⋘";
var llcorner = "⌞";
var Lleftarrow = "⇚";
var llhard = "⥫";
var lltri = "◺";
var Lmidot = "Ŀ";
var lmidot = "ŀ";
var lmoustache = "⎰";
var lmoust = "⎰";
var lnap = "⪉";
var lnapprox = "⪉";
var lne = "⪇";
var lnE = "≨";
var lneq = "⪇";
var lneqq = "≨";
var lnsim = "⋦";
var loang = "⟬";
var loarr = "⇽";
var lobrk = "⟦";
var longleftarrow = "⟵";
var LongLeftArrow = "⟵";
var Longleftarrow = "⟸";
var longleftrightarrow = "⟷";
var LongLeftRightArrow = "⟷";
var Longleftrightarrow = "⟺";
var longmapsto = "⟼";
var longrightarrow = "⟶";
var LongRightArrow = "⟶";
var Longrightarrow = "⟹";
var looparrowleft = "↫";
var looparrowright = "↬";
var lopar = "⦅";
var Lopf = "𝕃";
var lopf = "𝕝";
var loplus = "⨭";
var lotimes = "⨴";
var lowast = "∗";
var lowbar = "_";
var LowerLeftArrow = "↙";
var LowerRightArrow = "↘";
var loz = "◊";
var lozenge = "◊";
var lozf = "⧫";
var lpar = "(";
var lparlt = "⦓";
var lrarr = "⇆";
var lrcorner = "⌟";
var lrhar = "⇋";
var lrhard = "⥭";
var lrm = "‎";
var lrtri = "⊿";
var lsaquo = "‹";
var lscr = "𝓁";
var Lscr = "ℒ";
var lsh = "↰";
var Lsh = "↰";
var lsim = "≲";
var lsime = "⪍";
var lsimg = "⪏";
var lsqb = "[";
var lsquo = "‘";
var lsquor = "‚";
var Lstrok = "Ł";
var lstrok = "ł";
var ltcc = "⪦";
var ltcir = "⩹";
var lt$2 = "<";
var LT$1 = "<";
var Lt = "≪";
var ltdot = "⋖";
var lthree = "⋋";
var ltimes = "⋉";
var ltlarr = "⥶";
var ltquest = "⩻";
var ltri = "◃";
var ltrie = "⊴";
var ltrif = "◂";
var ltrPar = "⦖";
var lurdshar = "⥊";
var luruhar = "⥦";
var lvertneqq = "≨︀";
var lvnE = "≨︀";
var macr$1 = "¯";
var male = "♂";
var malt = "✠";
var maltese = "✠";
var map = "↦";
var mapsto = "↦";
var mapstodown = "↧";
var mapstoleft = "↤";
var mapstoup = "↥";
var marker = "▮";
var mcomma = "⨩";
var Mcy = "М";
var mcy = "м";
var mdash = "—";
var mDDot = "∺";
var measuredangle = "∡";
var MediumSpace = " ";
var Mellintrf = "ℳ";
var Mfr = "𝔐";
var mfr = "𝔪";
var mho = "℧";
var micro$1 = "µ";
var midast = "*";
var midcir = "⫰";
var mid = "∣";
var middot$1 = "·";
var minusb = "⊟";
var minus = "−";
var minusd = "∸";
var minusdu = "⨪";
var MinusPlus = "∓";
var mlcp = "⫛";
var mldr = "…";
var mnplus = "∓";
var models = "⊧";
var Mopf = "𝕄";
var mopf = "𝕞";
var mp = "∓";
var mscr = "𝓂";
var Mscr = "ℳ";
var mstpos = "∾";
var Mu = "Μ";
var mu = "μ";
var multimap = "⊸";
var mumap = "⊸";
var nabla = "∇";
var Nacute = "Ń";
var nacute = "ń";
var nang = "∠⃒";
var nap = "≉";
var napE = "⩰̸";
var napid = "≋̸";
var napos = "ŉ";
var napprox = "≉";
var natural = "♮";
var naturals = "ℕ";
var natur = "♮";
var nbsp$1 = " ";
var nbump = "≎̸";
var nbumpe = "≏̸";
var ncap = "⩃";
var Ncaron = "Ň";
var ncaron = "ň";
var Ncedil = "Ņ";
var ncedil = "ņ";
var ncong = "≇";
var ncongdot = "⩭̸";
var ncup = "⩂";
var Ncy = "Н";
var ncy = "н";
var ndash = "–";
var nearhk = "⤤";
var nearr = "↗";
var neArr = "⇗";
var nearrow = "↗";
var ne = "≠";
var nedot = "≐̸";
var NegativeMediumSpace = "​";
var NegativeThickSpace = "​";
var NegativeThinSpace = "​";
var NegativeVeryThinSpace = "​";
var nequiv = "≢";
var nesear = "⤨";
var nesim = "≂̸";
var NestedGreaterGreater = "≫";
var NestedLessLess = "≪";
var NewLine = "\n";
var nexist = "∄";
var nexists = "∄";
var Nfr = "𝔑";
var nfr = "𝔫";
var ngE = "≧̸";
var nge = "≱";
var ngeq = "≱";
var ngeqq = "≧̸";
var ngeqslant = "⩾̸";
var nges = "⩾̸";
var nGg = "⋙̸";
var ngsim = "≵";
var nGt = "≫⃒";
var ngt = "≯";
var ngtr = "≯";
var nGtv = "≫̸";
var nharr = "↮";
var nhArr = "⇎";
var nhpar = "⫲";
var ni = "∋";
var nis = "⋼";
var nisd = "⋺";
var niv = "∋";
var NJcy = "Њ";
var njcy = "њ";
var nlarr = "↚";
var nlArr = "⇍";
var nldr = "‥";
var nlE = "≦̸";
var nle = "≰";
var nleftarrow = "↚";
var nLeftarrow = "⇍";
var nleftrightarrow = "↮";
var nLeftrightarrow = "⇎";
var nleq = "≰";
var nleqq = "≦̸";
var nleqslant = "⩽̸";
var nles = "⩽̸";
var nless = "≮";
var nLl = "⋘̸";
var nlsim = "≴";
var nLt = "≪⃒";
var nlt = "≮";
var nltri = "⋪";
var nltrie = "⋬";
var nLtv = "≪̸";
var nmid = "∤";
var NoBreak = "⁠";
var NonBreakingSpace = " ";
var nopf = "𝕟";
var Nopf = "ℕ";
var Not = "⫬";
var not$1 = "¬";
var NotCongruent = "≢";
var NotCupCap = "≭";
var NotDoubleVerticalBar = "∦";
var NotElement = "∉";
var NotEqual = "≠";
var NotEqualTilde = "≂̸";
var NotExists = "∄";
var NotGreater = "≯";
var NotGreaterEqual = "≱";
var NotGreaterFullEqual = "≧̸";
var NotGreaterGreater = "≫̸";
var NotGreaterLess = "≹";
var NotGreaterSlantEqual = "⩾̸";
var NotGreaterTilde = "≵";
var NotHumpDownHump = "≎̸";
var NotHumpEqual = "≏̸";
var notin = "∉";
var notindot = "⋵̸";
var notinE = "⋹̸";
var notinva = "∉";
var notinvb = "⋷";
var notinvc = "⋶";
var NotLeftTriangleBar = "⧏̸";
var NotLeftTriangle = "⋪";
var NotLeftTriangleEqual = "⋬";
var NotLess = "≮";
var NotLessEqual = "≰";
var NotLessGreater = "≸";
var NotLessLess = "≪̸";
var NotLessSlantEqual = "⩽̸";
var NotLessTilde = "≴";
var NotNestedGreaterGreater = "⪢̸";
var NotNestedLessLess = "⪡̸";
var notni = "∌";
var notniva = "∌";
var notnivb = "⋾";
var notnivc = "⋽";
var NotPrecedes = "⊀";
var NotPrecedesEqual = "⪯̸";
var NotPrecedesSlantEqual = "⋠";
var NotReverseElement = "∌";
var NotRightTriangleBar = "⧐̸";
var NotRightTriangle = "⋫";
var NotRightTriangleEqual = "⋭";
var NotSquareSubset = "⊏̸";
var NotSquareSubsetEqual = "⋢";
var NotSquareSuperset = "⊐̸";
var NotSquareSupersetEqual = "⋣";
var NotSubset = "⊂⃒";
var NotSubsetEqual = "⊈";
var NotSucceeds = "⊁";
var NotSucceedsEqual = "⪰̸";
var NotSucceedsSlantEqual = "⋡";
var NotSucceedsTilde = "≿̸";
var NotSuperset = "⊃⃒";
var NotSupersetEqual = "⊉";
var NotTilde = "≁";
var NotTildeEqual = "≄";
var NotTildeFullEqual = "≇";
var NotTildeTilde = "≉";
var NotVerticalBar = "∤";
var nparallel = "∦";
var npar = "∦";
var nparsl = "⫽⃥";
var npart = "∂̸";
var npolint = "⨔";
var npr = "⊀";
var nprcue = "⋠";
var nprec = "⊀";
var npreceq = "⪯̸";
var npre = "⪯̸";
var nrarrc = "⤳̸";
var nrarr = "↛";
var nrArr = "⇏";
var nrarrw = "↝̸";
var nrightarrow = "↛";
var nRightarrow = "⇏";
var nrtri = "⋫";
var nrtrie = "⋭";
var nsc = "⊁";
var nsccue = "⋡";
var nsce = "⪰̸";
var Nscr = "𝒩";
var nscr = "𝓃";
var nshortmid = "∤";
var nshortparallel = "∦";
var nsim = "≁";
var nsime = "≄";
var nsimeq = "≄";
var nsmid = "∤";
var nspar = "∦";
var nsqsube = "⋢";
var nsqsupe = "⋣";
var nsub = "⊄";
var nsubE = "⫅̸";
var nsube = "⊈";
var nsubset = "⊂⃒";
var nsubseteq = "⊈";
var nsubseteqq = "⫅̸";
var nsucc = "⊁";
var nsucceq = "⪰̸";
var nsup = "⊅";
var nsupE = "⫆̸";
var nsupe = "⊉";
var nsupset = "⊃⃒";
var nsupseteq = "⊉";
var nsupseteqq = "⫆̸";
var ntgl = "≹";
var Ntilde$1 = "Ñ";
var ntilde$1 = "ñ";
var ntlg = "≸";
var ntriangleleft = "⋪";
var ntrianglelefteq = "⋬";
var ntriangleright = "⋫";
var ntrianglerighteq = "⋭";
var Nu = "Ν";
var nu = "ν";
var num = "#";
var numero = "№";
var numsp = " ";
var nvap = "≍⃒";
var nvdash = "⊬";
var nvDash = "⊭";
var nVdash = "⊮";
var nVDash = "⊯";
var nvge = "≥⃒";
var nvgt = ">⃒";
var nvHarr = "⤄";
var nvinfin = "⧞";
var nvlArr = "⤂";
var nvle = "≤⃒";
var nvlt = "<⃒";
var nvltrie = "⊴⃒";
var nvrArr = "⤃";
var nvrtrie = "⊵⃒";
var nvsim = "∼⃒";
var nwarhk = "⤣";
var nwarr = "↖";
var nwArr = "⇖";
var nwarrow = "↖";
var nwnear = "⤧";
var Oacute$1 = "Ó";
var oacute$1 = "ó";
var oast = "⊛";
var Ocirc$1 = "Ô";
var ocirc$1 = "ô";
var ocir = "⊚";
var Ocy = "О";
var ocy = "о";
var odash = "⊝";
var Odblac = "Ő";
var odblac = "ő";
var odiv = "⨸";
var odot = "⊙";
var odsold = "⦼";
var OElig = "Œ";
var oelig = "œ";
var ofcir = "⦿";
var Ofr = "𝔒";
var ofr = "𝔬";
var ogon = "˛";
var Ograve$1 = "Ò";
var ograve$1 = "ò";
var ogt = "⧁";
var ohbar = "⦵";
var ohm = "Ω";
var oint = "∮";
var olarr = "↺";
var olcir = "⦾";
var olcross = "⦻";
var oline = "‾";
var olt = "⧀";
var Omacr = "Ō";
var omacr = "ō";
var Omega = "Ω";
var omega = "ω";
var Omicron = "Ο";
var omicron = "ο";
var omid = "⦶";
var ominus = "⊖";
var Oopf = "𝕆";
var oopf = "𝕠";
var opar = "⦷";
var OpenCurlyDoubleQuote = "“";
var OpenCurlyQuote = "‘";
var operp = "⦹";
var oplus = "⊕";
var orarr = "↻";
var Or = "⩔";
var or = "∨";
var ord = "⩝";
var order = "ℴ";
var orderof = "ℴ";
var ordf$1 = "ª";
var ordm$1 = "º";
var origof = "⊶";
var oror = "⩖";
var orslope = "⩗";
var orv = "⩛";
var oS = "Ⓢ";
var Oscr = "𝒪";
var oscr = "ℴ";
var Oslash$1 = "Ø";
var oslash$1 = "ø";
var osol = "⊘";
var Otilde$1 = "Õ";
var otilde$1 = "õ";
var otimesas = "⨶";
var Otimes = "⨷";
var otimes = "⊗";
var Ouml$1 = "Ö";
var ouml$1 = "ö";
var ovbar = "⌽";
var OverBar = "‾";
var OverBrace = "⏞";
var OverBracket = "⎴";
var OverParenthesis = "⏜";
var para$1 = "¶";
var parallel = "∥";
var par = "∥";
var parsim = "⫳";
var parsl = "⫽";
var part = "∂";
var PartialD = "∂";
var Pcy = "П";
var pcy = "п";
var percnt = "%";
var period = ".";
var permil = "‰";
var perp = "⊥";
var pertenk = "‱";
var Pfr = "𝔓";
var pfr = "𝔭";
var Phi = "Φ";
var phi = "φ";
var phiv = "ϕ";
var phmmat = "ℳ";
var phone = "☎";
var Pi = "Π";
var pi = "π";
var pitchfork = "⋔";
var piv = "ϖ";
var planck = "ℏ";
var planckh = "ℎ";
var plankv = "ℏ";
var plusacir = "⨣";
var plusb = "⊞";
var pluscir = "⨢";
var plus$1 = "+";
var plusdo = "∔";
var plusdu = "⨥";
var pluse = "⩲";
var PlusMinus = "±";
var plusmn$1 = "±";
var plussim = "⨦";
var plustwo = "⨧";
var pm = "±";
var Poincareplane = "ℌ";
var pointint = "⨕";
var popf = "𝕡";
var Popf = "ℙ";
var pound$1 = "£";
var prap = "⪷";
var Pr = "⪻";
var pr = "≺";
var prcue = "≼";
var precapprox = "⪷";
var prec = "≺";
var preccurlyeq = "≼";
var Precedes = "≺";
var PrecedesEqual = "⪯";
var PrecedesSlantEqual = "≼";
var PrecedesTilde = "≾";
var preceq = "⪯";
var precnapprox = "⪹";
var precneqq = "⪵";
var precnsim = "⋨";
var pre = "⪯";
var prE = "⪳";
var precsim = "≾";
var prime = "′";
var Prime = "″";
var primes = "ℙ";
var prnap = "⪹";
var prnE = "⪵";
var prnsim = "⋨";
var prod = "∏";
var Product = "∏";
var profalar = "⌮";
var profline = "⌒";
var profsurf = "⌓";
var prop = "∝";
var Proportional = "∝";
var Proportion = "∷";
var propto = "∝";
var prsim = "≾";
var prurel = "⊰";
var Pscr = "𝒫";
var pscr = "𝓅";
var Psi = "Ψ";
var psi = "ψ";
var puncsp = " ";
var Qfr = "𝔔";
var qfr = "𝔮";
var qint = "⨌";
var qopf = "𝕢";
var Qopf = "ℚ";
var qprime = "⁗";
var Qscr = "𝒬";
var qscr = "𝓆";
var quaternions = "ℍ";
var quatint = "⨖";
var quest = "?";
var questeq = "≟";
var quot$2 = "\"";
var QUOT$1 = "\"";
var rAarr = "⇛";
var race = "∽̱";
var Racute = "Ŕ";
var racute = "ŕ";
var radic = "√";
var raemptyv = "⦳";
var rang = "⟩";
var Rang = "⟫";
var rangd = "⦒";
var range = "⦥";
var rangle = "⟩";
var raquo$1 = "»";
var rarrap = "⥵";
var rarrb = "⇥";
var rarrbfs = "⤠";
var rarrc = "⤳";
var rarr = "→";
var Rarr = "↠";
var rArr = "⇒";
var rarrfs = "⤞";
var rarrhk = "↪";
var rarrlp = "↬";
var rarrpl = "⥅";
var rarrsim = "⥴";
var Rarrtl = "⤖";
var rarrtl = "↣";
var rarrw = "↝";
var ratail = "⤚";
var rAtail = "⤜";
var ratio = "∶";
var rationals = "ℚ";
var rbarr = "⤍";
var rBarr = "⤏";
var RBarr = "⤐";
var rbbrk = "❳";
var rbrace = "}";
var rbrack = "]";
var rbrke = "⦌";
var rbrksld = "⦎";
var rbrkslu = "⦐";
var Rcaron = "Ř";
var rcaron = "ř";
var Rcedil = "Ŗ";
var rcedil = "ŗ";
var rceil = "⌉";
var rcub = "}";
var Rcy = "Р";
var rcy = "р";
var rdca = "⤷";
var rdldhar = "⥩";
var rdquo = "”";
var rdquor = "”";
var rdsh = "↳";
var real = "ℜ";
var realine = "ℛ";
var realpart = "ℜ";
var reals = "ℝ";
var Re = "ℜ";
var rect = "▭";
var reg$1 = "®";
var REG$1 = "®";
var ReverseElement = "∋";
var ReverseEquilibrium = "⇋";
var ReverseUpEquilibrium = "⥯";
var rfisht = "⥽";
var rfloor = "⌋";
var rfr = "𝔯";
var Rfr = "ℜ";
var rHar = "⥤";
var rhard = "⇁";
var rharu = "⇀";
var rharul = "⥬";
var Rho = "Ρ";
var rho = "ρ";
var rhov = "ϱ";
var RightAngleBracket = "⟩";
var RightArrowBar = "⇥";
var rightarrow = "→";
var RightArrow = "→";
var Rightarrow = "⇒";
var RightArrowLeftArrow = "⇄";
var rightarrowtail = "↣";
var RightCeiling = "⌉";
var RightDoubleBracket = "⟧";
var RightDownTeeVector = "⥝";
var RightDownVectorBar = "⥕";
var RightDownVector = "⇂";
var RightFloor = "⌋";
var rightharpoondown = "⇁";
var rightharpoonup = "⇀";
var rightleftarrows = "⇄";
var rightleftharpoons = "⇌";
var rightrightarrows = "⇉";
var rightsquigarrow = "↝";
var RightTeeArrow = "↦";
var RightTee = "⊢";
var RightTeeVector = "⥛";
var rightthreetimes = "⋌";
var RightTriangleBar = "⧐";
var RightTriangle = "⊳";
var RightTriangleEqual = "⊵";
var RightUpDownVector = "⥏";
var RightUpTeeVector = "⥜";
var RightUpVectorBar = "⥔";
var RightUpVector = "↾";
var RightVectorBar = "⥓";
var RightVector = "⇀";
var ring = "˚";
var risingdotseq = "≓";
var rlarr = "⇄";
var rlhar = "⇌";
var rlm = "‏";
var rmoustache = "⎱";
var rmoust = "⎱";
var rnmid = "⫮";
var roang = "⟭";
var roarr = "⇾";
var robrk = "⟧";
var ropar = "⦆";
var ropf = "𝕣";
var Ropf = "ℝ";
var roplus = "⨮";
var rotimes = "⨵";
var RoundImplies = "⥰";
var rpar = ")";
var rpargt = "⦔";
var rppolint = "⨒";
var rrarr = "⇉";
var Rrightarrow = "⇛";
var rsaquo = "›";
var rscr = "𝓇";
var Rscr = "ℛ";
var rsh = "↱";
var Rsh = "↱";
var rsqb = "]";
var rsquo = "’";
var rsquor = "’";
var rthree = "⋌";
var rtimes = "⋊";
var rtri = "▹";
var rtrie = "⊵";
var rtrif = "▸";
var rtriltri = "⧎";
var RuleDelayed = "⧴";
var ruluhar = "⥨";
var rx = "℞";
var Sacute = "Ś";
var sacute = "ś";
var sbquo = "‚";
var scap = "⪸";
var Scaron = "Š";
var scaron = "š";
var Sc = "⪼";
var sc = "≻";
var sccue = "≽";
var sce = "⪰";
var scE = "⪴";
var Scedil = "Ş";
var scedil = "ş";
var Scirc = "Ŝ";
var scirc = "ŝ";
var scnap = "⪺";
var scnE = "⪶";
var scnsim = "⋩";
var scpolint = "⨓";
var scsim = "≿";
var Scy = "С";
var scy = "с";
var sdotb = "⊡";
var sdot = "⋅";
var sdote = "⩦";
var searhk = "⤥";
var searr = "↘";
var seArr = "⇘";
var searrow = "↘";
var sect$1 = "§";
var semi = ";";
var seswar = "⤩";
var setminus = "∖";
var setmn = "∖";
var sext = "✶";
var Sfr = "𝔖";
var sfr = "𝔰";
var sfrown = "⌢";
var sharp = "♯";
var SHCHcy = "Щ";
var shchcy = "щ";
var SHcy = "Ш";
var shcy = "ш";
var ShortDownArrow = "↓";
var ShortLeftArrow = "←";
var shortmid = "∣";
var shortparallel = "∥";
var ShortRightArrow = "→";
var ShortUpArrow = "↑";
var shy$1 = "­";
var Sigma = "Σ";
var sigma = "σ";
var sigmaf = "ς";
var sigmav = "ς";
var sim = "∼";
var simdot = "⩪";
var sime = "≃";
var simeq = "≃";
var simg = "⪞";
var simgE = "⪠";
var siml = "⪝";
var simlE = "⪟";
var simne = "≆";
var simplus = "⨤";
var simrarr = "⥲";
var slarr = "←";
var SmallCircle = "∘";
var smallsetminus = "∖";
var smashp = "⨳";
var smeparsl = "⧤";
var smid = "∣";
var smile = "⌣";
var smt = "⪪";
var smte = "⪬";
var smtes = "⪬︀";
var SOFTcy = "Ь";
var softcy = "ь";
var solbar = "⌿";
var solb = "⧄";
var sol = "/";
var Sopf = "𝕊";
var sopf = "𝕤";
var spades = "♠";
var spadesuit = "♠";
var spar = "∥";
var sqcap = "⊓";
var sqcaps = "⊓︀";
var sqcup = "⊔";
var sqcups = "⊔︀";
var Sqrt = "√";
var sqsub = "⊏";
var sqsube = "⊑";
var sqsubset = "⊏";
var sqsubseteq = "⊑";
var sqsup = "⊐";
var sqsupe = "⊒";
var sqsupset = "⊐";
var sqsupseteq = "⊒";
var square = "□";
var Square = "□";
var SquareIntersection = "⊓";
var SquareSubset = "⊏";
var SquareSubsetEqual = "⊑";
var SquareSuperset = "⊐";
var SquareSupersetEqual = "⊒";
var SquareUnion = "⊔";
var squarf = "▪";
var squ = "□";
var squf = "▪";
var srarr = "→";
var Sscr = "𝒮";
var sscr = "𝓈";
var ssetmn = "∖";
var ssmile = "⌣";
var sstarf = "⋆";
var Star = "⋆";
var star = "☆";
var starf = "★";
var straightepsilon = "ϵ";
var straightphi = "ϕ";
var strns = "¯";
var sub = "⊂";
var Sub = "⋐";
var subdot = "⪽";
var subE = "⫅";
var sube = "⊆";
var subedot = "⫃";
var submult = "⫁";
var subnE = "⫋";
var subne = "⊊";
var subplus = "⪿";
var subrarr = "⥹";
var subset = "⊂";
var Subset = "⋐";
var subseteq = "⊆";
var subseteqq = "⫅";
var SubsetEqual = "⊆";
var subsetneq = "⊊";
var subsetneqq = "⫋";
var subsim = "⫇";
var subsub = "⫕";
var subsup = "⫓";
var succapprox = "⪸";
var succ = "≻";
var succcurlyeq = "≽";
var Succeeds = "≻";
var SucceedsEqual = "⪰";
var SucceedsSlantEqual = "≽";
var SucceedsTilde = "≿";
var succeq = "⪰";
var succnapprox = "⪺";
var succneqq = "⪶";
var succnsim = "⋩";
var succsim = "≿";
var SuchThat = "∋";
var sum = "∑";
var Sum = "∑";
var sung = "♪";
var sup1$1 = "¹";
var sup2$1 = "²";
var sup3$1 = "³";
var sup = "⊃";
var Sup = "⋑";
var supdot = "⪾";
var supdsub = "⫘";
var supE = "⫆";
var supe = "⊇";
var supedot = "⫄";
var Superset = "⊃";
var SupersetEqual = "⊇";
var suphsol = "⟉";
var suphsub = "⫗";
var suplarr = "⥻";
var supmult = "⫂";
var supnE = "⫌";
var supne = "⊋";
var supplus = "⫀";
var supset = "⊃";
var Supset = "⋑";
var supseteq = "⊇";
var supseteqq = "⫆";
var supsetneq = "⊋";
var supsetneqq = "⫌";
var supsim = "⫈";
var supsub = "⫔";
var supsup = "⫖";
var swarhk = "⤦";
var swarr = "↙";
var swArr = "⇙";
var swarrow = "↙";
var swnwar = "⤪";
var szlig$1 = "ß";
var Tab = "\t";
var target = "⌖";
var Tau = "Τ";
var tau = "τ";
var tbrk = "⎴";
var Tcaron = "Ť";
var tcaron = "ť";
var Tcedil = "Ţ";
var tcedil = "ţ";
var Tcy = "Т";
var tcy = "т";
var tdot = "⃛";
var telrec = "⌕";
var Tfr = "𝔗";
var tfr = "𝔱";
var there4 = "∴";
var therefore = "∴";
var Therefore = "∴";
var Theta = "Θ";
var theta = "θ";
var thetasym = "ϑ";
var thetav = "ϑ";
var thickapprox = "≈";
var thicksim = "∼";
var ThickSpace = "  ";
var ThinSpace = " ";
var thinsp = " ";
var thkap = "≈";
var thksim = "∼";
var THORN$1 = "Þ";
var thorn$1 = "þ";
var tilde$1 = "˜";
var Tilde = "∼";
var TildeEqual = "≃";
var TildeFullEqual = "≅";
var TildeTilde = "≈";
var timesbar = "⨱";
var timesb = "⊠";
var times$1 = "×";
var timesd = "⨰";
var tint = "∭";
var toea = "⤨";
var topbot = "⌶";
var topcir = "⫱";
var top = "⊤";
var Topf = "𝕋";
var topf = "𝕥";
var topfork = "⫚";
var tosa = "⤩";
var tprime = "‴";
var trade = "™";
var TRADE = "™";
var triangle = "▵";
var triangledown = "▿";
var triangleleft = "◃";
var trianglelefteq = "⊴";
var triangleq = "≜";
var triangleright = "▹";
var trianglerighteq = "⊵";
var tridot = "◬";
var trie = "≜";
var triminus = "⨺";
var TripleDot = "⃛";
var triplus = "⨹";
var trisb = "⧍";
var tritime = "⨻";
var trpezium = "⏢";
var Tscr = "𝒯";
var tscr = "𝓉";
var TScy = "Ц";
var tscy = "ц";
var TSHcy = "Ћ";
var tshcy = "ћ";
var Tstrok = "Ŧ";
var tstrok = "ŧ";
var twixt = "≬";
var twoheadleftarrow = "↞";
var twoheadrightarrow = "↠";
var Uacute$1 = "Ú";
var uacute$1 = "ú";
var uarr = "↑";
var Uarr = "↟";
var uArr = "⇑";
var Uarrocir = "⥉";
var Ubrcy = "Ў";
var ubrcy = "ў";
var Ubreve = "Ŭ";
var ubreve = "ŭ";
var Ucirc$1 = "Û";
var ucirc$1 = "û";
var Ucy = "У";
var ucy = "у";
var udarr = "⇅";
var Udblac = "Ű";
var udblac = "ű";
var udhar = "⥮";
var ufisht = "⥾";
var Ufr = "𝔘";
var ufr = "𝔲";
var Ugrave$1 = "Ù";
var ugrave$1 = "ù";
var uHar = "⥣";
var uharl = "↿";
var uharr = "↾";
var uhblk = "▀";
var ulcorn = "⌜";
var ulcorner = "⌜";
var ulcrop = "⌏";
var ultri = "◸";
var Umacr = "Ū";
var umacr = "ū";
var uml$1 = "¨";
var UnderBar = "_";
var UnderBrace = "⏟";
var UnderBracket = "⎵";
var UnderParenthesis = "⏝";
var Union = "⋃";
var UnionPlus = "⊎";
var Uogon = "Ų";
var uogon = "ų";
var Uopf = "𝕌";
var uopf = "𝕦";
var UpArrowBar = "⤒";
var uparrow = "↑";
var UpArrow = "↑";
var Uparrow = "⇑";
var UpArrowDownArrow = "⇅";
var updownarrow = "↕";
var UpDownArrow = "↕";
var Updownarrow = "⇕";
var UpEquilibrium = "⥮";
var upharpoonleft = "↿";
var upharpoonright = "↾";
var uplus = "⊎";
var UpperLeftArrow = "↖";
var UpperRightArrow = "↗";
var upsi = "υ";
var Upsi = "ϒ";
var upsih = "ϒ";
var Upsilon = "Υ";
var upsilon = "υ";
var UpTeeArrow = "↥";
var UpTee = "⊥";
var upuparrows = "⇈";
var urcorn = "⌝";
var urcorner = "⌝";
var urcrop = "⌎";
var Uring = "Ů";
var uring = "ů";
var urtri = "◹";
var Uscr = "𝒰";
var uscr = "𝓊";
var utdot = "⋰";
var Utilde = "Ũ";
var utilde = "ũ";
var utri = "▵";
var utrif = "▴";
var uuarr = "⇈";
var Uuml$1 = "Ü";
var uuml$1 = "ü";
var uwangle = "⦧";
var vangrt = "⦜";
var varepsilon = "ϵ";
var varkappa = "ϰ";
var varnothing = "∅";
var varphi = "ϕ";
var varpi = "ϖ";
var varpropto = "∝";
var varr = "↕";
var vArr = "⇕";
var varrho = "ϱ";
var varsigma = "ς";
var varsubsetneq = "⊊︀";
var varsubsetneqq = "⫋︀";
var varsupsetneq = "⊋︀";
var varsupsetneqq = "⫌︀";
var vartheta = "ϑ";
var vartriangleleft = "⊲";
var vartriangleright = "⊳";
var vBar = "⫨";
var Vbar = "⫫";
var vBarv = "⫩";
var Vcy = "В";
var vcy = "в";
var vdash = "⊢";
var vDash = "⊨";
var Vdash = "⊩";
var VDash = "⊫";
var Vdashl = "⫦";
var veebar = "⊻";
var vee = "∨";
var Vee = "⋁";
var veeeq = "≚";
var vellip = "⋮";
var verbar = "|";
var Verbar = "‖";
var vert = "|";
var Vert = "‖";
var VerticalBar = "∣";
var VerticalLine = "|";
var VerticalSeparator = "❘";
var VerticalTilde = "≀";
var VeryThinSpace = " ";
var Vfr = "𝔙";
var vfr = "𝔳";
var vltri = "⊲";
var vnsub = "⊂⃒";
var vnsup = "⊃⃒";
var Vopf = "𝕍";
var vopf = "𝕧";
var vprop = "∝";
var vrtri = "⊳";
var Vscr = "𝒱";
var vscr = "𝓋";
var vsubnE = "⫋︀";
var vsubne = "⊊︀";
var vsupnE = "⫌︀";
var vsupne = "⊋︀";
var Vvdash = "⊪";
var vzigzag = "⦚";
var Wcirc = "Ŵ";
var wcirc = "ŵ";
var wedbar = "⩟";
var wedge = "∧";
var Wedge = "⋀";
var wedgeq = "≙";
var weierp = "℘";
var Wfr = "𝔚";
var wfr = "𝔴";
var Wopf = "𝕎";
var wopf = "𝕨";
var wp = "℘";
var wr = "≀";
var wreath = "≀";
var Wscr = "𝒲";
var wscr = "𝓌";
var xcap = "⋂";
var xcirc = "◯";
var xcup = "⋃";
var xdtri = "▽";
var Xfr = "𝔛";
var xfr = "𝔵";
var xharr = "⟷";
var xhArr = "⟺";
var Xi = "Ξ";
var xi = "ξ";
var xlarr = "⟵";
var xlArr = "⟸";
var xmap = "⟼";
var xnis = "⋻";
var xodot = "⨀";
var Xopf = "𝕏";
var xopf = "𝕩";
var xoplus = "⨁";
var xotime = "⨂";
var xrarr = "⟶";
var xrArr = "⟹";
var Xscr = "𝒳";
var xscr = "𝓍";
var xsqcup = "⨆";
var xuplus = "⨄";
var xutri = "△";
var xvee = "⋁";
var xwedge = "⋀";
var Yacute$1 = "Ý";
var yacute$1 = "ý";
var YAcy = "Я";
var yacy = "я";
var Ycirc = "Ŷ";
var ycirc = "ŷ";
var Ycy = "Ы";
var ycy = "ы";
var yen$1 = "¥";
var Yfr = "𝔜";
var yfr = "𝔶";
var YIcy = "Ї";
var yicy = "ї";
var Yopf = "𝕐";
var yopf = "𝕪";
var Yscr = "𝒴";
var yscr = "𝓎";
var YUcy = "Ю";
var yucy = "ю";
var yuml$1 = "ÿ";
var Yuml = "Ÿ";
var Zacute = "Ź";
var zacute = "ź";
var Zcaron = "Ž";
var zcaron = "ž";
var Zcy = "З";
var zcy = "з";
var Zdot = "Ż";
var zdot = "ż";
var zeetrf = "ℨ";
var ZeroWidthSpace = "​";
var Zeta = "Ζ";
var zeta = "ζ";
var zfr = "𝔷";
var Zfr = "ℨ";
var ZHcy = "Ж";
var zhcy = "ж";
var zigrarr = "⇝";
var zopf = "𝕫";
var Zopf = "ℤ";
var Zscr = "𝒵";
var zscr = "𝓏";
var zwj = "‍";
var zwnj = "‌";
var require$$1$1 = {
	Aacute: Aacute$1,
	aacute: aacute$1,
	Abreve: Abreve,
	abreve: abreve,
	ac: ac,
	acd: acd,
	acE: acE,
	Acirc: Acirc$1,
	acirc: acirc$1,
	acute: acute$1,
	Acy: Acy,
	acy: acy,
	AElig: AElig$1,
	aelig: aelig$1,
	af: af,
	Afr: Afr,
	afr: afr,
	Agrave: Agrave$1,
	agrave: agrave$1,
	alefsym: alefsym,
	aleph: aleph,
	Alpha: Alpha,
	alpha: alpha,
	Amacr: Amacr,
	amacr: amacr,
	amalg: amalg,
	amp: amp$2,
	AMP: AMP$1,
	andand: andand,
	And: And,
	and: and,
	andd: andd,
	andslope: andslope,
	andv: andv,
	ang: ang,
	ange: ange,
	angle: angle,
	angmsdaa: angmsdaa,
	angmsdab: angmsdab,
	angmsdac: angmsdac,
	angmsdad: angmsdad,
	angmsdae: angmsdae,
	angmsdaf: angmsdaf,
	angmsdag: angmsdag,
	angmsdah: angmsdah,
	angmsd: angmsd,
	angrt: angrt,
	angrtvb: angrtvb,
	angrtvbd: angrtvbd,
	angsph: angsph,
	angst: angst,
	angzarr: angzarr,
	Aogon: Aogon,
	aogon: aogon,
	Aopf: Aopf,
	aopf: aopf,
	apacir: apacir,
	ap: ap,
	apE: apE,
	ape: ape,
	apid: apid,
	apos: apos$1,
	ApplyFunction: ApplyFunction,
	approx: approx,
	approxeq: approxeq,
	Aring: Aring$1,
	aring: aring$1,
	Ascr: Ascr,
	ascr: ascr,
	Assign: Assign,
	ast: ast,
	asymp: asymp,
	asympeq: asympeq,
	Atilde: Atilde$1,
	atilde: atilde$1,
	Auml: Auml$1,
	auml: auml$1,
	awconint: awconint,
	awint: awint,
	backcong: backcong,
	backepsilon: backepsilon,
	backprime: backprime,
	backsim: backsim,
	backsimeq: backsimeq,
	Backslash: Backslash,
	Barv: Barv,
	barvee: barvee,
	barwed: barwed,
	Barwed: Barwed,
	barwedge: barwedge,
	bbrk: bbrk,
	bbrktbrk: bbrktbrk,
	bcong: bcong,
	Bcy: Bcy,
	bcy: bcy,
	bdquo: bdquo,
	becaus: becaus,
	because: because,
	Because: Because,
	bemptyv: bemptyv,
	bepsi: bepsi,
	bernou: bernou,
	Bernoullis: Bernoullis,
	Beta: Beta,
	beta: beta,
	beth: beth,
	between: between,
	Bfr: Bfr,
	bfr: bfr,
	bigcap: bigcap,
	bigcirc: bigcirc,
	bigcup: bigcup,
	bigodot: bigodot,
	bigoplus: bigoplus,
	bigotimes: bigotimes,
	bigsqcup: bigsqcup,
	bigstar: bigstar,
	bigtriangledown: bigtriangledown,
	bigtriangleup: bigtriangleup,
	biguplus: biguplus,
	bigvee: bigvee,
	bigwedge: bigwedge,
	bkarow: bkarow,
	blacklozenge: blacklozenge,
	blacksquare: blacksquare,
	blacktriangle: blacktriangle,
	blacktriangledown: blacktriangledown,
	blacktriangleleft: blacktriangleleft,
	blacktriangleright: blacktriangleright,
	blank: blank,
	blk12: blk12,
	blk14: blk14,
	blk34: blk34,
	block: block,
	bne: bne,
	bnequiv: bnequiv,
	bNot: bNot,
	bnot: bnot,
	Bopf: Bopf,
	bopf: bopf,
	bot: bot,
	bottom: bottom,
	bowtie: bowtie,
	boxbox: boxbox,
	boxdl: boxdl,
	boxdL: boxdL,
	boxDl: boxDl,
	boxDL: boxDL,
	boxdr: boxdr,
	boxdR: boxdR,
	boxDr: boxDr,
	boxDR: boxDR,
	boxh: boxh,
	boxH: boxH,
	boxhd: boxhd,
	boxHd: boxHd,
	boxhD: boxhD,
	boxHD: boxHD,
	boxhu: boxhu,
	boxHu: boxHu,
	boxhU: boxhU,
	boxHU: boxHU,
	boxminus: boxminus,
	boxplus: boxplus,
	boxtimes: boxtimes,
	boxul: boxul,
	boxuL: boxuL,
	boxUl: boxUl,
	boxUL: boxUL,
	boxur: boxur,
	boxuR: boxuR,
	boxUr: boxUr,
	boxUR: boxUR,
	boxv: boxv,
	boxV: boxV,
	boxvh: boxvh,
	boxvH: boxvH,
	boxVh: boxVh,
	boxVH: boxVH,
	boxvl: boxvl,
	boxvL: boxvL,
	boxVl: boxVl,
	boxVL: boxVL,
	boxvr: boxvr,
	boxvR: boxvR,
	boxVr: boxVr,
	boxVR: boxVR,
	bprime: bprime,
	breve: breve,
	Breve: Breve,
	brvbar: brvbar$1,
	bscr: bscr,
	Bscr: Bscr,
	bsemi: bsemi,
	bsim: bsim,
	bsime: bsime,
	bsolb: bsolb,
	bsol: bsol,
	bsolhsub: bsolhsub,
	bull: bull,
	bullet: bullet,
	bump: bump,
	bumpE: bumpE,
	bumpe: bumpe,
	Bumpeq: Bumpeq,
	bumpeq: bumpeq,
	Cacute: Cacute,
	cacute: cacute,
	capand: capand,
	capbrcup: capbrcup,
	capcap: capcap,
	cap: cap,
	Cap: Cap,
	capcup: capcup,
	capdot: capdot,
	CapitalDifferentialD: CapitalDifferentialD,
	caps: caps,
	caret: caret$1,
	caron: caron,
	Cayleys: Cayleys,
	ccaps: ccaps,
	Ccaron: Ccaron,
	ccaron: ccaron,
	Ccedil: Ccedil$1,
	ccedil: ccedil$1,
	Ccirc: Ccirc,
	ccirc: ccirc,
	Cconint: Cconint,
	ccups: ccups,
	ccupssm: ccupssm,
	Cdot: Cdot,
	cdot: cdot,
	cedil: cedil$1,
	Cedilla: Cedilla,
	cemptyv: cemptyv,
	cent: cent$1,
	centerdot: centerdot,
	CenterDot: CenterDot,
	cfr: cfr,
	Cfr: Cfr,
	CHcy: CHcy,
	chcy: chcy,
	check: check,
	checkmark: checkmark,
	Chi: Chi,
	chi: chi,
	circ: circ,
	circeq: circeq,
	circlearrowleft: circlearrowleft,
	circlearrowright: circlearrowright,
	circledast: circledast,
	circledcirc: circledcirc,
	circleddash: circleddash,
	CircleDot: CircleDot,
	circledR: circledR,
	circledS: circledS,
	CircleMinus: CircleMinus,
	CirclePlus: CirclePlus,
	CircleTimes: CircleTimes,
	cir: cir,
	cirE: cirE,
	cire: cire,
	cirfnint: cirfnint,
	cirmid: cirmid,
	cirscir: cirscir,
	ClockwiseContourIntegral: ClockwiseContourIntegral,
	CloseCurlyDoubleQuote: CloseCurlyDoubleQuote,
	CloseCurlyQuote: CloseCurlyQuote,
	clubs: clubs,
	clubsuit: clubsuit,
	colon: colon,
	Colon: Colon,
	Colone: Colone,
	colone: colone,
	coloneq: coloneq,
	comma: comma$1,
	commat: commat,
	comp: comp,
	compfn: compfn,
	complement: complement,
	complexes: complexes,
	cong: cong,
	congdot: congdot,
	Congruent: Congruent,
	conint: conint,
	Conint: Conint,
	ContourIntegral: ContourIntegral,
	copf: copf,
	Copf: Copf,
	coprod: coprod,
	Coproduct: Coproduct,
	copy: copy$1,
	COPY: COPY$1,
	copysr: copysr,
	CounterClockwiseContourIntegral: CounterClockwiseContourIntegral,
	crarr: crarr,
	cross: cross,
	Cross: Cross,
	Cscr: Cscr,
	cscr: cscr,
	csub: csub,
	csube: csube,
	csup: csup,
	csupe: csupe,
	ctdot: ctdot,
	cudarrl: cudarrl,
	cudarrr: cudarrr,
	cuepr: cuepr,
	cuesc: cuesc,
	cularr: cularr,
	cularrp: cularrp,
	cupbrcap: cupbrcap,
	cupcap: cupcap,
	CupCap: CupCap,
	cup: cup,
	Cup: Cup,
	cupcup: cupcup,
	cupdot: cupdot,
	cupor: cupor,
	cups: cups,
	curarr: curarr,
	curarrm: curarrm,
	curlyeqprec: curlyeqprec,
	curlyeqsucc: curlyeqsucc,
	curlyvee: curlyvee,
	curlywedge: curlywedge,
	curren: curren$1,
	curvearrowleft: curvearrowleft,
	curvearrowright: curvearrowright,
	cuvee: cuvee,
	cuwed: cuwed,
	cwconint: cwconint,
	cwint: cwint,
	cylcty: cylcty,
	dagger: dagger,
	Dagger: Dagger,
	daleth: daleth,
	darr: darr,
	Darr: Darr,
	dArr: dArr,
	dash: dash,
	Dashv: Dashv,
	dashv: dashv,
	dbkarow: dbkarow,
	dblac: dblac,
	Dcaron: Dcaron,
	dcaron: dcaron,
	Dcy: Dcy,
	dcy: dcy,
	ddagger: ddagger,
	ddarr: ddarr,
	DD: DD,
	dd: dd,
	DDotrahd: DDotrahd,
	ddotseq: ddotseq,
	deg: deg$1,
	Del: Del,
	Delta: Delta,
	delta: delta,
	demptyv: demptyv,
	dfisht: dfisht,
	Dfr: Dfr,
	dfr: dfr,
	dHar: dHar,
	dharl: dharl,
	dharr: dharr,
	DiacriticalAcute: DiacriticalAcute,
	DiacriticalDot: DiacriticalDot,
	DiacriticalDoubleAcute: DiacriticalDoubleAcute,
	DiacriticalGrave: DiacriticalGrave,
	DiacriticalTilde: DiacriticalTilde,
	diam: diam,
	diamond: diamond,
	Diamond: Diamond,
	diamondsuit: diamondsuit,
	diams: diams,
	die: die,
	DifferentialD: DifferentialD,
	digamma: digamma,
	disin: disin,
	div: div,
	divide: divide$1,
	divideontimes: divideontimes,
	divonx: divonx,
	DJcy: DJcy,
	djcy: djcy,
	dlcorn: dlcorn,
	dlcrop: dlcrop,
	dollar: dollar$1,
	Dopf: Dopf,
	dopf: dopf,
	Dot: Dot,
	dot: dot,
	DotDot: DotDot,
	doteq: doteq,
	doteqdot: doteqdot,
	DotEqual: DotEqual,
	dotminus: dotminus,
	dotplus: dotplus,
	dotsquare: dotsquare,
	doublebarwedge: doublebarwedge,
	DoubleContourIntegral: DoubleContourIntegral,
	DoubleDot: DoubleDot,
	DoubleDownArrow: DoubleDownArrow,
	DoubleLeftArrow: DoubleLeftArrow,
	DoubleLeftRightArrow: DoubleLeftRightArrow,
	DoubleLeftTee: DoubleLeftTee,
	DoubleLongLeftArrow: DoubleLongLeftArrow,
	DoubleLongLeftRightArrow: DoubleLongLeftRightArrow,
	DoubleLongRightArrow: DoubleLongRightArrow,
	DoubleRightArrow: DoubleRightArrow,
	DoubleRightTee: DoubleRightTee,
	DoubleUpArrow: DoubleUpArrow,
	DoubleUpDownArrow: DoubleUpDownArrow,
	DoubleVerticalBar: DoubleVerticalBar,
	DownArrowBar: DownArrowBar,
	downarrow: downarrow,
	DownArrow: DownArrow,
	Downarrow: Downarrow,
	DownArrowUpArrow: DownArrowUpArrow,
	DownBreve: DownBreve,
	downdownarrows: downdownarrows,
	downharpoonleft: downharpoonleft,
	downharpoonright: downharpoonright,
	DownLeftRightVector: DownLeftRightVector,
	DownLeftTeeVector: DownLeftTeeVector,
	DownLeftVectorBar: DownLeftVectorBar,
	DownLeftVector: DownLeftVector,
	DownRightTeeVector: DownRightTeeVector,
	DownRightVectorBar: DownRightVectorBar,
	DownRightVector: DownRightVector,
	DownTeeArrow: DownTeeArrow,
	DownTee: DownTee,
	drbkarow: drbkarow,
	drcorn: drcorn,
	drcrop: drcrop,
	Dscr: Dscr,
	dscr: dscr,
	DScy: DScy,
	dscy: dscy,
	dsol: dsol,
	Dstrok: Dstrok,
	dstrok: dstrok,
	dtdot: dtdot,
	dtri: dtri,
	dtrif: dtrif,
	duarr: duarr,
	duhar: duhar,
	dwangle: dwangle,
	DZcy: DZcy,
	dzcy: dzcy,
	dzigrarr: dzigrarr,
	Eacute: Eacute$1,
	eacute: eacute$1,
	easter: easter,
	Ecaron: Ecaron,
	ecaron: ecaron,
	Ecirc: Ecirc$1,
	ecirc: ecirc$1,
	ecir: ecir,
	ecolon: ecolon,
	Ecy: Ecy,
	ecy: ecy,
	eDDot: eDDot,
	Edot: Edot,
	edot: edot,
	eDot: eDot,
	ee: ee,
	efDot: efDot,
	Efr: Efr,
	efr: efr,
	eg: eg,
	Egrave: Egrave$1,
	egrave: egrave$1,
	egs: egs,
	egsdot: egsdot,
	el: el,
	Element: Element,
	elinters: elinters,
	ell: ell,
	els: els,
	elsdot: elsdot,
	Emacr: Emacr,
	emacr: emacr,
	empty: empty,
	emptyset: emptyset,
	EmptySmallSquare: EmptySmallSquare,
	emptyv: emptyv,
	EmptyVerySmallSquare: EmptyVerySmallSquare,
	emsp13: emsp13,
	emsp14: emsp14,
	emsp: emsp,
	ENG: ENG,
	eng: eng,
	ensp: ensp,
	Eogon: Eogon,
	eogon: eogon,
	Eopf: Eopf,
	eopf: eopf,
	epar: epar,
	eparsl: eparsl,
	eplus: eplus,
	epsi: epsi,
	Epsilon: Epsilon,
	epsilon: epsilon,
	epsiv: epsiv,
	eqcirc: eqcirc,
	eqcolon: eqcolon,
	eqsim: eqsim,
	eqslantgtr: eqslantgtr,
	eqslantless: eqslantless,
	Equal: Equal,
	equals: equals,
	EqualTilde: EqualTilde,
	equest: equest,
	Equilibrium: Equilibrium,
	equiv: equiv,
	equivDD: equivDD,
	eqvparsl: eqvparsl,
	erarr: erarr,
	erDot: erDot,
	escr: escr,
	Escr: Escr,
	esdot: esdot,
	Esim: Esim,
	esim: esim,
	Eta: Eta,
	eta: eta,
	ETH: ETH$1,
	eth: eth$1,
	Euml: Euml$1,
	euml: euml$1,
	euro: euro,
	excl: excl,
	exist: exist,
	Exists: Exists,
	expectation: expectation,
	exponentiale: exponentiale,
	ExponentialE: ExponentialE,
	fallingdotseq: fallingdotseq,
	Fcy: Fcy,
	fcy: fcy,
	female: female,
	ffilig: ffilig,
	fflig: fflig,
	ffllig: ffllig,
	Ffr: Ffr,
	ffr: ffr,
	filig: filig,
	FilledSmallSquare: FilledSmallSquare,
	FilledVerySmallSquare: FilledVerySmallSquare,
	fjlig: fjlig,
	flat: flat,
	fllig: fllig,
	fltns: fltns,
	fnof: fnof,
	Fopf: Fopf,
	fopf: fopf,
	forall: forall,
	ForAll: ForAll,
	fork: fork,
	forkv: forkv,
	Fouriertrf: Fouriertrf,
	fpartint: fpartint,
	frac12: frac12$1,
	frac13: frac13,
	frac14: frac14$1,
	frac15: frac15,
	frac16: frac16,
	frac18: frac18,
	frac23: frac23,
	frac25: frac25,
	frac34: frac34$1,
	frac35: frac35,
	frac38: frac38,
	frac45: frac45,
	frac56: frac56,
	frac58: frac58,
	frac78: frac78,
	frasl: frasl,
	frown: frown,
	fscr: fscr,
	Fscr: Fscr,
	gacute: gacute,
	Gamma: Gamma,
	gamma: gamma,
	Gammad: Gammad,
	gammad: gammad,
	gap: gap,
	Gbreve: Gbreve,
	gbreve: gbreve,
	Gcedil: Gcedil,
	Gcirc: Gcirc,
	gcirc: gcirc,
	Gcy: Gcy,
	gcy: gcy,
	Gdot: Gdot,
	gdot: gdot,
	ge: ge,
	gE: gE,
	gEl: gEl,
	gel: gel,
	geq: geq,
	geqq: geqq,
	geqslant: geqslant,
	gescc: gescc,
	ges: ges,
	gesdot: gesdot,
	gesdoto: gesdoto,
	gesdotol: gesdotol,
	gesl: gesl,
	gesles: gesles,
	Gfr: Gfr,
	gfr: gfr,
	gg: gg,
	Gg: Gg,
	ggg: ggg,
	gimel: gimel,
	GJcy: GJcy,
	gjcy: gjcy,
	gla: gla,
	gl: gl,
	glE: glE,
	glj: glj,
	gnap: gnap,
	gnapprox: gnapprox,
	gne: gne,
	gnE: gnE,
	gneq: gneq,
	gneqq: gneqq,
	gnsim: gnsim,
	Gopf: Gopf,
	gopf: gopf,
	grave: grave,
	GreaterEqual: GreaterEqual,
	GreaterEqualLess: GreaterEqualLess,
	GreaterFullEqual: GreaterFullEqual,
	GreaterGreater: GreaterGreater,
	GreaterLess: GreaterLess,
	GreaterSlantEqual: GreaterSlantEqual,
	GreaterTilde: GreaterTilde,
	Gscr: Gscr,
	gscr: gscr,
	gsim: gsim,
	gsime: gsime,
	gsiml: gsiml,
	gtcc: gtcc,
	gtcir: gtcir,
	gt: gt$3,
	GT: GT$1,
	Gt: Gt,
	gtdot: gtdot,
	gtlPar: gtlPar,
	gtquest: gtquest,
	gtrapprox: gtrapprox,
	gtrarr: gtrarr,
	gtrdot: gtrdot,
	gtreqless: gtreqless,
	gtreqqless: gtreqqless,
	gtrless: gtrless,
	gtrsim: gtrsim,
	gvertneqq: gvertneqq,
	gvnE: gvnE,
	Hacek: Hacek,
	hairsp: hairsp,
	half: half,
	hamilt: hamilt,
	HARDcy: HARDcy,
	hardcy: hardcy,
	harrcir: harrcir,
	harr: harr,
	hArr: hArr,
	harrw: harrw,
	Hat: Hat,
	hbar: hbar,
	Hcirc: Hcirc,
	hcirc: hcirc,
	hearts: hearts,
	heartsuit: heartsuit,
	hellip: hellip,
	hercon: hercon,
	hfr: hfr,
	Hfr: Hfr,
	HilbertSpace: HilbertSpace,
	hksearow: hksearow,
	hkswarow: hkswarow,
	hoarr: hoarr,
	homtht: homtht,
	hookleftarrow: hookleftarrow,
	hookrightarrow: hookrightarrow,
	hopf: hopf,
	Hopf: Hopf,
	horbar: horbar,
	HorizontalLine: HorizontalLine,
	hscr: hscr,
	Hscr: Hscr,
	hslash: hslash,
	Hstrok: Hstrok,
	hstrok: hstrok,
	HumpDownHump: HumpDownHump,
	HumpEqual: HumpEqual,
	hybull: hybull,
	hyphen: hyphen,
	Iacute: Iacute$1,
	iacute: iacute$1,
	ic: ic,
	Icirc: Icirc$1,
	icirc: icirc$1,
	Icy: Icy,
	icy: icy,
	Idot: Idot,
	IEcy: IEcy,
	iecy: iecy,
	iexcl: iexcl$1,
	iff: iff,
	ifr: ifr,
	Ifr: Ifr,
	Igrave: Igrave$1,
	igrave: igrave$1,
	ii: ii,
	iiiint: iiiint,
	iiint: iiint,
	iinfin: iinfin,
	iiota: iiota,
	IJlig: IJlig,
	ijlig: ijlig,
	Imacr: Imacr,
	imacr: imacr,
	image: image,
	ImaginaryI: ImaginaryI,
	imagline: imagline,
	imagpart: imagpart,
	imath: imath,
	Im: Im,
	imof: imof,
	imped: imped,
	Implies: Implies,
	incare: incare,
	"in": "∈",
	infin: infin,
	infintie: infintie,
	inodot: inodot,
	intcal: intcal,
	int: int,
	Int: Int,
	integers: integers,
	Integral: Integral,
	intercal: intercal,
	Intersection: Intersection,
	intlarhk: intlarhk,
	intprod: intprod,
	InvisibleComma: InvisibleComma,
	InvisibleTimes: InvisibleTimes,
	IOcy: IOcy,
	iocy: iocy,
	Iogon: Iogon,
	iogon: iogon,
	Iopf: Iopf,
	iopf: iopf,
	Iota: Iota,
	iota: iota,
	iprod: iprod,
	iquest: iquest$1,
	iscr: iscr,
	Iscr: Iscr,
	isin: isin,
	isindot: isindot,
	isinE: isinE,
	isins: isins,
	isinsv: isinsv,
	isinv: isinv,
	it: it,
	Itilde: Itilde,
	itilde: itilde,
	Iukcy: Iukcy,
	iukcy: iukcy,
	Iuml: Iuml$1,
	iuml: iuml$1,
	Jcirc: Jcirc,
	jcirc: jcirc,
	Jcy: Jcy,
	jcy: jcy,
	Jfr: Jfr,
	jfr: jfr,
	jmath: jmath,
	Jopf: Jopf,
	jopf: jopf,
	Jscr: Jscr,
	jscr: jscr,
	Jsercy: Jsercy,
	jsercy: jsercy,
	Jukcy: Jukcy,
	jukcy: jukcy,
	Kappa: Kappa,
	kappa: kappa,
	kappav: kappav,
	Kcedil: Kcedil,
	kcedil: kcedil,
	Kcy: Kcy,
	kcy: kcy,
	Kfr: Kfr,
	kfr: kfr,
	kgreen: kgreen,
	KHcy: KHcy,
	khcy: khcy,
	KJcy: KJcy,
	kjcy: kjcy,
	Kopf: Kopf,
	kopf: kopf,
	Kscr: Kscr,
	kscr: kscr,
	lAarr: lAarr,
	Lacute: Lacute,
	lacute: lacute,
	laemptyv: laemptyv,
	lagran: lagran,
	Lambda: Lambda,
	lambda: lambda,
	lang: lang,
	Lang: Lang,
	langd: langd,
	langle: langle,
	lap: lap,
	Laplacetrf: Laplacetrf,
	laquo: laquo$1,
	larrb: larrb,
	larrbfs: larrbfs,
	larr: larr,
	Larr: Larr,
	lArr: lArr,
	larrfs: larrfs,
	larrhk: larrhk,
	larrlp: larrlp,
	larrpl: larrpl,
	larrsim: larrsim,
	larrtl: larrtl,
	latail: latail,
	lAtail: lAtail,
	lat: lat,
	late: late,
	lates: lates,
	lbarr: lbarr,
	lBarr: lBarr,
	lbbrk: lbbrk,
	lbrace: lbrace,
	lbrack: lbrack,
	lbrke: lbrke,
	lbrksld: lbrksld,
	lbrkslu: lbrkslu,
	Lcaron: Lcaron,
	lcaron: lcaron,
	Lcedil: Lcedil,
	lcedil: lcedil,
	lceil: lceil,
	lcub: lcub,
	Lcy: Lcy,
	lcy: lcy,
	ldca: ldca,
	ldquo: ldquo,
	ldquor: ldquor,
	ldrdhar: ldrdhar,
	ldrushar: ldrushar,
	ldsh: ldsh,
	le: le,
	lE: lE,
	LeftAngleBracket: LeftAngleBracket,
	LeftArrowBar: LeftArrowBar,
	leftarrow: leftarrow,
	LeftArrow: LeftArrow,
	Leftarrow: Leftarrow,
	LeftArrowRightArrow: LeftArrowRightArrow,
	leftarrowtail: leftarrowtail,
	LeftCeiling: LeftCeiling,
	LeftDoubleBracket: LeftDoubleBracket,
	LeftDownTeeVector: LeftDownTeeVector,
	LeftDownVectorBar: LeftDownVectorBar,
	LeftDownVector: LeftDownVector,
	LeftFloor: LeftFloor,
	leftharpoondown: leftharpoondown,
	leftharpoonup: leftharpoonup,
	leftleftarrows: leftleftarrows,
	leftrightarrow: leftrightarrow,
	LeftRightArrow: LeftRightArrow,
	Leftrightarrow: Leftrightarrow,
	leftrightarrows: leftrightarrows,
	leftrightharpoons: leftrightharpoons,
	leftrightsquigarrow: leftrightsquigarrow,
	LeftRightVector: LeftRightVector,
	LeftTeeArrow: LeftTeeArrow,
	LeftTee: LeftTee,
	LeftTeeVector: LeftTeeVector,
	leftthreetimes: leftthreetimes,
	LeftTriangleBar: LeftTriangleBar,
	LeftTriangle: LeftTriangle,
	LeftTriangleEqual: LeftTriangleEqual,
	LeftUpDownVector: LeftUpDownVector,
	LeftUpTeeVector: LeftUpTeeVector,
	LeftUpVectorBar: LeftUpVectorBar,
	LeftUpVector: LeftUpVector,
	LeftVectorBar: LeftVectorBar,
	LeftVector: LeftVector,
	lEg: lEg,
	leg: leg,
	leq: leq,
	leqq: leqq,
	leqslant: leqslant,
	lescc: lescc,
	les: les,
	lesdot: lesdot,
	lesdoto: lesdoto,
	lesdotor: lesdotor,
	lesg: lesg,
	lesges: lesges,
	lessapprox: lessapprox,
	lessdot: lessdot,
	lesseqgtr: lesseqgtr,
	lesseqqgtr: lesseqqgtr,
	LessEqualGreater: LessEqualGreater,
	LessFullEqual: LessFullEqual,
	LessGreater: LessGreater,
	lessgtr: lessgtr,
	LessLess: LessLess,
	lesssim: lesssim,
	LessSlantEqual: LessSlantEqual,
	LessTilde: LessTilde,
	lfisht: lfisht,
	lfloor: lfloor,
	Lfr: Lfr,
	lfr: lfr,
	lg: lg,
	lgE: lgE,
	lHar: lHar,
	lhard: lhard,
	lharu: lharu,
	lharul: lharul,
	lhblk: lhblk,
	LJcy: LJcy,
	ljcy: ljcy,
	llarr: llarr,
	ll: ll,
	Ll: Ll,
	llcorner: llcorner,
	Lleftarrow: Lleftarrow,
	llhard: llhard,
	lltri: lltri,
	Lmidot: Lmidot,
	lmidot: lmidot,
	lmoustache: lmoustache,
	lmoust: lmoust,
	lnap: lnap,
	lnapprox: lnapprox,
	lne: lne,
	lnE: lnE,
	lneq: lneq,
	lneqq: lneqq,
	lnsim: lnsim,
	loang: loang,
	loarr: loarr,
	lobrk: lobrk,
	longleftarrow: longleftarrow,
	LongLeftArrow: LongLeftArrow,
	Longleftarrow: Longleftarrow,
	longleftrightarrow: longleftrightarrow,
	LongLeftRightArrow: LongLeftRightArrow,
	Longleftrightarrow: Longleftrightarrow,
	longmapsto: longmapsto,
	longrightarrow: longrightarrow,
	LongRightArrow: LongRightArrow,
	Longrightarrow: Longrightarrow,
	looparrowleft: looparrowleft,
	looparrowright: looparrowright,
	lopar: lopar,
	Lopf: Lopf,
	lopf: lopf,
	loplus: loplus,
	lotimes: lotimes,
	lowast: lowast,
	lowbar: lowbar,
	LowerLeftArrow: LowerLeftArrow,
	LowerRightArrow: LowerRightArrow,
	loz: loz,
	lozenge: lozenge,
	lozf: lozf,
	lpar: lpar,
	lparlt: lparlt,
	lrarr: lrarr,
	lrcorner: lrcorner,
	lrhar: lrhar,
	lrhard: lrhard,
	lrm: lrm,
	lrtri: lrtri,
	lsaquo: lsaquo,
	lscr: lscr,
	Lscr: Lscr,
	lsh: lsh,
	Lsh: Lsh,
	lsim: lsim,
	lsime: lsime,
	lsimg: lsimg,
	lsqb: lsqb,
	lsquo: lsquo,
	lsquor: lsquor,
	Lstrok: Lstrok,
	lstrok: lstrok,
	ltcc: ltcc,
	ltcir: ltcir,
	lt: lt$2,
	LT: LT$1,
	Lt: Lt,
	ltdot: ltdot,
	lthree: lthree,
	ltimes: ltimes,
	ltlarr: ltlarr,
	ltquest: ltquest,
	ltri: ltri,
	ltrie: ltrie,
	ltrif: ltrif,
	ltrPar: ltrPar,
	lurdshar: lurdshar,
	luruhar: luruhar,
	lvertneqq: lvertneqq,
	lvnE: lvnE,
	macr: macr$1,
	male: male,
	malt: malt,
	maltese: maltese,
	"Map": "⤅",
	map: map,
	mapsto: mapsto,
	mapstodown: mapstodown,
	mapstoleft: mapstoleft,
	mapstoup: mapstoup,
	marker: marker,
	mcomma: mcomma,
	Mcy: Mcy,
	mcy: mcy,
	mdash: mdash,
	mDDot: mDDot,
	measuredangle: measuredangle,
	MediumSpace: MediumSpace,
	Mellintrf: Mellintrf,
	Mfr: Mfr,
	mfr: mfr,
	mho: mho,
	micro: micro$1,
	midast: midast,
	midcir: midcir,
	mid: mid,
	middot: middot$1,
	minusb: minusb,
	minus: minus,
	minusd: minusd,
	minusdu: minusdu,
	MinusPlus: MinusPlus,
	mlcp: mlcp,
	mldr: mldr,
	mnplus: mnplus,
	models: models,
	Mopf: Mopf,
	mopf: mopf,
	mp: mp,
	mscr: mscr,
	Mscr: Mscr,
	mstpos: mstpos,
	Mu: Mu,
	mu: mu,
	multimap: multimap,
	mumap: mumap,
	nabla: nabla,
	Nacute: Nacute,
	nacute: nacute,
	nang: nang,
	nap: nap,
	napE: napE,
	napid: napid,
	napos: napos,
	napprox: napprox,
	natural: natural,
	naturals: naturals,
	natur: natur,
	nbsp: nbsp$1,
	nbump: nbump,
	nbumpe: nbumpe,
	ncap: ncap,
	Ncaron: Ncaron,
	ncaron: ncaron,
	Ncedil: Ncedil,
	ncedil: ncedil,
	ncong: ncong,
	ncongdot: ncongdot,
	ncup: ncup,
	Ncy: Ncy,
	ncy: ncy,
	ndash: ndash,
	nearhk: nearhk,
	nearr: nearr,
	neArr: neArr,
	nearrow: nearrow,
	ne: ne,
	nedot: nedot,
	NegativeMediumSpace: NegativeMediumSpace,
	NegativeThickSpace: NegativeThickSpace,
	NegativeThinSpace: NegativeThinSpace,
	NegativeVeryThinSpace: NegativeVeryThinSpace,
	nequiv: nequiv,
	nesear: nesear,
	nesim: nesim,
	NestedGreaterGreater: NestedGreaterGreater,
	NestedLessLess: NestedLessLess,
	NewLine: NewLine,
	nexist: nexist,
	nexists: nexists,
	Nfr: Nfr,
	nfr: nfr,
	ngE: ngE,
	nge: nge,
	ngeq: ngeq,
	ngeqq: ngeqq,
	ngeqslant: ngeqslant,
	nges: nges,
	nGg: nGg,
	ngsim: ngsim,
	nGt: nGt,
	ngt: ngt,
	ngtr: ngtr,
	nGtv: nGtv,
	nharr: nharr,
	nhArr: nhArr,
	nhpar: nhpar,
	ni: ni,
	nis: nis,
	nisd: nisd,
	niv: niv,
	NJcy: NJcy,
	njcy: njcy,
	nlarr: nlarr,
	nlArr: nlArr,
	nldr: nldr,
	nlE: nlE,
	nle: nle,
	nleftarrow: nleftarrow,
	nLeftarrow: nLeftarrow,
	nleftrightarrow: nleftrightarrow,
	nLeftrightarrow: nLeftrightarrow,
	nleq: nleq,
	nleqq: nleqq,
	nleqslant: nleqslant,
	nles: nles,
	nless: nless,
	nLl: nLl,
	nlsim: nlsim,
	nLt: nLt,
	nlt: nlt,
	nltri: nltri,
	nltrie: nltrie,
	nLtv: nLtv,
	nmid: nmid,
	NoBreak: NoBreak,
	NonBreakingSpace: NonBreakingSpace,
	nopf: nopf,
	Nopf: Nopf,
	Not: Not,
	not: not$1,
	NotCongruent: NotCongruent,
	NotCupCap: NotCupCap,
	NotDoubleVerticalBar: NotDoubleVerticalBar,
	NotElement: NotElement,
	NotEqual: NotEqual,
	NotEqualTilde: NotEqualTilde,
	NotExists: NotExists,
	NotGreater: NotGreater,
	NotGreaterEqual: NotGreaterEqual,
	NotGreaterFullEqual: NotGreaterFullEqual,
	NotGreaterGreater: NotGreaterGreater,
	NotGreaterLess: NotGreaterLess,
	NotGreaterSlantEqual: NotGreaterSlantEqual,
	NotGreaterTilde: NotGreaterTilde,
	NotHumpDownHump: NotHumpDownHump,
	NotHumpEqual: NotHumpEqual,
	notin: notin,
	notindot: notindot,
	notinE: notinE,
	notinva: notinva,
	notinvb: notinvb,
	notinvc: notinvc,
	NotLeftTriangleBar: NotLeftTriangleBar,
	NotLeftTriangle: NotLeftTriangle,
	NotLeftTriangleEqual: NotLeftTriangleEqual,
	NotLess: NotLess,
	NotLessEqual: NotLessEqual,
	NotLessGreater: NotLessGreater,
	NotLessLess: NotLessLess,
	NotLessSlantEqual: NotLessSlantEqual,
	NotLessTilde: NotLessTilde,
	NotNestedGreaterGreater: NotNestedGreaterGreater,
	NotNestedLessLess: NotNestedLessLess,
	notni: notni,
	notniva: notniva,
	notnivb: notnivb,
	notnivc: notnivc,
	NotPrecedes: NotPrecedes,
	NotPrecedesEqual: NotPrecedesEqual,
	NotPrecedesSlantEqual: NotPrecedesSlantEqual,
	NotReverseElement: NotReverseElement,
	NotRightTriangleBar: NotRightTriangleBar,
	NotRightTriangle: NotRightTriangle,
	NotRightTriangleEqual: NotRightTriangleEqual,
	NotSquareSubset: NotSquareSubset,
	NotSquareSubsetEqual: NotSquareSubsetEqual,
	NotSquareSuperset: NotSquareSuperset,
	NotSquareSupersetEqual: NotSquareSupersetEqual,
	NotSubset: NotSubset,
	NotSubsetEqual: NotSubsetEqual,
	NotSucceeds: NotSucceeds,
	NotSucceedsEqual: NotSucceedsEqual,
	NotSucceedsSlantEqual: NotSucceedsSlantEqual,
	NotSucceedsTilde: NotSucceedsTilde,
	NotSuperset: NotSuperset,
	NotSupersetEqual: NotSupersetEqual,
	NotTilde: NotTilde,
	NotTildeEqual: NotTildeEqual,
	NotTildeFullEqual: NotTildeFullEqual,
	NotTildeTilde: NotTildeTilde,
	NotVerticalBar: NotVerticalBar,
	nparallel: nparallel,
	npar: npar,
	nparsl: nparsl,
	npart: npart,
	npolint: npolint,
	npr: npr,
	nprcue: nprcue,
	nprec: nprec,
	npreceq: npreceq,
	npre: npre,
	nrarrc: nrarrc,
	nrarr: nrarr,
	nrArr: nrArr,
	nrarrw: nrarrw,
	nrightarrow: nrightarrow,
	nRightarrow: nRightarrow,
	nrtri: nrtri,
	nrtrie: nrtrie,
	nsc: nsc,
	nsccue: nsccue,
	nsce: nsce,
	Nscr: Nscr,
	nscr: nscr,
	nshortmid: nshortmid,
	nshortparallel: nshortparallel,
	nsim: nsim,
	nsime: nsime,
	nsimeq: nsimeq,
	nsmid: nsmid,
	nspar: nspar,
	nsqsube: nsqsube,
	nsqsupe: nsqsupe,
	nsub: nsub,
	nsubE: nsubE,
	nsube: nsube,
	nsubset: nsubset,
	nsubseteq: nsubseteq,
	nsubseteqq: nsubseteqq,
	nsucc: nsucc,
	nsucceq: nsucceq,
	nsup: nsup,
	nsupE: nsupE,
	nsupe: nsupe,
	nsupset: nsupset,
	nsupseteq: nsupseteq,
	nsupseteqq: nsupseteqq,
	ntgl: ntgl,
	Ntilde: Ntilde$1,
	ntilde: ntilde$1,
	ntlg: ntlg,
	ntriangleleft: ntriangleleft,
	ntrianglelefteq: ntrianglelefteq,
	ntriangleright: ntriangleright,
	ntrianglerighteq: ntrianglerighteq,
	Nu: Nu,
	nu: nu,
	num: num,
	numero: numero,
	numsp: numsp,
	nvap: nvap,
	nvdash: nvdash,
	nvDash: nvDash,
	nVdash: nVdash,
	nVDash: nVDash,
	nvge: nvge,
	nvgt: nvgt,
	nvHarr: nvHarr,
	nvinfin: nvinfin,
	nvlArr: nvlArr,
	nvle: nvle,
	nvlt: nvlt,
	nvltrie: nvltrie,
	nvrArr: nvrArr,
	nvrtrie: nvrtrie,
	nvsim: nvsim,
	nwarhk: nwarhk,
	nwarr: nwarr,
	nwArr: nwArr,
	nwarrow: nwarrow,
	nwnear: nwnear,
	Oacute: Oacute$1,
	oacute: oacute$1,
	oast: oast,
	Ocirc: Ocirc$1,
	ocirc: ocirc$1,
	ocir: ocir,
	Ocy: Ocy,
	ocy: ocy,
	odash: odash,
	Odblac: Odblac,
	odblac: odblac,
	odiv: odiv,
	odot: odot,
	odsold: odsold,
	OElig: OElig,
	oelig: oelig,
	ofcir: ofcir,
	Ofr: Ofr,
	ofr: ofr,
	ogon: ogon,
	Ograve: Ograve$1,
	ograve: ograve$1,
	ogt: ogt,
	ohbar: ohbar,
	ohm: ohm,
	oint: oint,
	olarr: olarr,
	olcir: olcir,
	olcross: olcross,
	oline: oline,
	olt: olt,
	Omacr: Omacr,
	omacr: omacr,
	Omega: Omega,
	omega: omega,
	Omicron: Omicron,
	omicron: omicron,
	omid: omid,
	ominus: ominus,
	Oopf: Oopf,
	oopf: oopf,
	opar: opar,
	OpenCurlyDoubleQuote: OpenCurlyDoubleQuote,
	OpenCurlyQuote: OpenCurlyQuote,
	operp: operp,
	oplus: oplus,
	orarr: orarr,
	Or: Or,
	or: or,
	ord: ord,
	order: order,
	orderof: orderof,
	ordf: ordf$1,
	ordm: ordm$1,
	origof: origof,
	oror: oror,
	orslope: orslope,
	orv: orv,
	oS: oS,
	Oscr: Oscr,
	oscr: oscr,
	Oslash: Oslash$1,
	oslash: oslash$1,
	osol: osol,
	Otilde: Otilde$1,
	otilde: otilde$1,
	otimesas: otimesas,
	Otimes: Otimes,
	otimes: otimes,
	Ouml: Ouml$1,
	ouml: ouml$1,
	ovbar: ovbar,
	OverBar: OverBar,
	OverBrace: OverBrace,
	OverBracket: OverBracket,
	OverParenthesis: OverParenthesis,
	para: para$1,
	parallel: parallel,
	par: par,
	parsim: parsim,
	parsl: parsl,
	part: part,
	PartialD: PartialD,
	Pcy: Pcy,
	pcy: pcy,
	percnt: percnt,
	period: period,
	permil: permil,
	perp: perp,
	pertenk: pertenk,
	Pfr: Pfr,
	pfr: pfr,
	Phi: Phi,
	phi: phi,
	phiv: phiv,
	phmmat: phmmat,
	phone: phone,
	Pi: Pi,
	pi: pi,
	pitchfork: pitchfork,
	piv: piv,
	planck: planck,
	planckh: planckh,
	plankv: plankv,
	plusacir: plusacir,
	plusb: plusb,
	pluscir: pluscir,
	plus: plus$1,
	plusdo: plusdo,
	plusdu: plusdu,
	pluse: pluse,
	PlusMinus: PlusMinus,
	plusmn: plusmn$1,
	plussim: plussim,
	plustwo: plustwo,
	pm: pm,
	Poincareplane: Poincareplane,
	pointint: pointint,
	popf: popf,
	Popf: Popf,
	pound: pound$1,
	prap: prap,
	Pr: Pr,
	pr: pr,
	prcue: prcue,
	precapprox: precapprox,
	prec: prec,
	preccurlyeq: preccurlyeq,
	Precedes: Precedes,
	PrecedesEqual: PrecedesEqual,
	PrecedesSlantEqual: PrecedesSlantEqual,
	PrecedesTilde: PrecedesTilde,
	preceq: preceq,
	precnapprox: precnapprox,
	precneqq: precneqq,
	precnsim: precnsim,
	pre: pre,
	prE: prE,
	precsim: precsim,
	prime: prime,
	Prime: Prime,
	primes: primes,
	prnap: prnap,
	prnE: prnE,
	prnsim: prnsim,
	prod: prod,
	Product: Product,
	profalar: profalar,
	profline: profline,
	profsurf: profsurf,
	prop: prop,
	Proportional: Proportional,
	Proportion: Proportion,
	propto: propto,
	prsim: prsim,
	prurel: prurel,
	Pscr: Pscr,
	pscr: pscr,
	Psi: Psi,
	psi: psi,
	puncsp: puncsp,
	Qfr: Qfr,
	qfr: qfr,
	qint: qint,
	qopf: qopf,
	Qopf: Qopf,
	qprime: qprime,
	Qscr: Qscr,
	qscr: qscr,
	quaternions: quaternions,
	quatint: quatint,
	quest: quest,
	questeq: questeq,
	quot: quot$2,
	QUOT: QUOT$1,
	rAarr: rAarr,
	race: race,
	Racute: Racute,
	racute: racute,
	radic: radic,
	raemptyv: raemptyv,
	rang: rang,
	Rang: Rang,
	rangd: rangd,
	range: range,
	rangle: rangle,
	raquo: raquo$1,
	rarrap: rarrap,
	rarrb: rarrb,
	rarrbfs: rarrbfs,
	rarrc: rarrc,
	rarr: rarr,
	Rarr: Rarr,
	rArr: rArr,
	rarrfs: rarrfs,
	rarrhk: rarrhk,
	rarrlp: rarrlp,
	rarrpl: rarrpl,
	rarrsim: rarrsim,
	Rarrtl: Rarrtl,
	rarrtl: rarrtl,
	rarrw: rarrw,
	ratail: ratail,
	rAtail: rAtail,
	ratio: ratio,
	rationals: rationals,
	rbarr: rbarr,
	rBarr: rBarr,
	RBarr: RBarr,
	rbbrk: rbbrk,
	rbrace: rbrace,
	rbrack: rbrack,
	rbrke: rbrke,
	rbrksld: rbrksld,
	rbrkslu: rbrkslu,
	Rcaron: Rcaron,
	rcaron: rcaron,
	Rcedil: Rcedil,
	rcedil: rcedil,
	rceil: rceil,
	rcub: rcub,
	Rcy: Rcy,
	rcy: rcy,
	rdca: rdca,
	rdldhar: rdldhar,
	rdquo: rdquo,
	rdquor: rdquor,
	rdsh: rdsh,
	real: real,
	realine: realine,
	realpart: realpart,
	reals: reals,
	Re: Re,
	rect: rect,
	reg: reg$1,
	REG: REG$1,
	ReverseElement: ReverseElement,
	ReverseEquilibrium: ReverseEquilibrium,
	ReverseUpEquilibrium: ReverseUpEquilibrium,
	rfisht: rfisht,
	rfloor: rfloor,
	rfr: rfr,
	Rfr: Rfr,
	rHar: rHar,
	rhard: rhard,
	rharu: rharu,
	rharul: rharul,
	Rho: Rho,
	rho: rho,
	rhov: rhov,
	RightAngleBracket: RightAngleBracket,
	RightArrowBar: RightArrowBar,
	rightarrow: rightarrow,
	RightArrow: RightArrow,
	Rightarrow: Rightarrow,
	RightArrowLeftArrow: RightArrowLeftArrow,
	rightarrowtail: rightarrowtail,
	RightCeiling: RightCeiling,
	RightDoubleBracket: RightDoubleBracket,
	RightDownTeeVector: RightDownTeeVector,
	RightDownVectorBar: RightDownVectorBar,
	RightDownVector: RightDownVector,
	RightFloor: RightFloor,
	rightharpoondown: rightharpoondown,
	rightharpoonup: rightharpoonup,
	rightleftarrows: rightleftarrows,
	rightleftharpoons: rightleftharpoons,
	rightrightarrows: rightrightarrows,
	rightsquigarrow: rightsquigarrow,
	RightTeeArrow: RightTeeArrow,
	RightTee: RightTee,
	RightTeeVector: RightTeeVector,
	rightthreetimes: rightthreetimes,
	RightTriangleBar: RightTriangleBar,
	RightTriangle: RightTriangle,
	RightTriangleEqual: RightTriangleEqual,
	RightUpDownVector: RightUpDownVector,
	RightUpTeeVector: RightUpTeeVector,
	RightUpVectorBar: RightUpVectorBar,
	RightUpVector: RightUpVector,
	RightVectorBar: RightVectorBar,
	RightVector: RightVector,
	ring: ring,
	risingdotseq: risingdotseq,
	rlarr: rlarr,
	rlhar: rlhar,
	rlm: rlm,
	rmoustache: rmoustache,
	rmoust: rmoust,
	rnmid: rnmid,
	roang: roang,
	roarr: roarr,
	robrk: robrk,
	ropar: ropar,
	ropf: ropf,
	Ropf: Ropf,
	roplus: roplus,
	rotimes: rotimes,
	RoundImplies: RoundImplies,
	rpar: rpar,
	rpargt: rpargt,
	rppolint: rppolint,
	rrarr: rrarr,
	Rrightarrow: Rrightarrow,
	rsaquo: rsaquo,
	rscr: rscr,
	Rscr: Rscr,
	rsh: rsh,
	Rsh: Rsh,
	rsqb: rsqb,
	rsquo: rsquo,
	rsquor: rsquor,
	rthree: rthree,
	rtimes: rtimes,
	rtri: rtri,
	rtrie: rtrie,
	rtrif: rtrif,
	rtriltri: rtriltri,
	RuleDelayed: RuleDelayed,
	ruluhar: ruluhar,
	rx: rx,
	Sacute: Sacute,
	sacute: sacute,
	sbquo: sbquo,
	scap: scap,
	Scaron: Scaron,
	scaron: scaron,
	Sc: Sc,
	sc: sc,
	sccue: sccue,
	sce: sce,
	scE: scE,
	Scedil: Scedil,
	scedil: scedil,
	Scirc: Scirc,
	scirc: scirc,
	scnap: scnap,
	scnE: scnE,
	scnsim: scnsim,
	scpolint: scpolint,
	scsim: scsim,
	Scy: Scy,
	scy: scy,
	sdotb: sdotb,
	sdot: sdot,
	sdote: sdote,
	searhk: searhk,
	searr: searr,
	seArr: seArr,
	searrow: searrow,
	sect: sect$1,
	semi: semi,
	seswar: seswar,
	setminus: setminus,
	setmn: setmn,
	sext: sext,
	Sfr: Sfr,
	sfr: sfr,
	sfrown: sfrown,
	sharp: sharp,
	SHCHcy: SHCHcy,
	shchcy: shchcy,
	SHcy: SHcy,
	shcy: shcy,
	ShortDownArrow: ShortDownArrow,
	ShortLeftArrow: ShortLeftArrow,
	shortmid: shortmid,
	shortparallel: shortparallel,
	ShortRightArrow: ShortRightArrow,
	ShortUpArrow: ShortUpArrow,
	shy: shy$1,
	Sigma: Sigma,
	sigma: sigma,
	sigmaf: sigmaf,
	sigmav: sigmav,
	sim: sim,
	simdot: simdot,
	sime: sime,
	simeq: simeq,
	simg: simg,
	simgE: simgE,
	siml: siml,
	simlE: simlE,
	simne: simne,
	simplus: simplus,
	simrarr: simrarr,
	slarr: slarr,
	SmallCircle: SmallCircle,
	smallsetminus: smallsetminus,
	smashp: smashp,
	smeparsl: smeparsl,
	smid: smid,
	smile: smile,
	smt: smt,
	smte: smte,
	smtes: smtes,
	SOFTcy: SOFTcy,
	softcy: softcy,
	solbar: solbar,
	solb: solb,
	sol: sol,
	Sopf: Sopf,
	sopf: sopf,
	spades: spades,
	spadesuit: spadesuit,
	spar: spar,
	sqcap: sqcap,
	sqcaps: sqcaps,
	sqcup: sqcup,
	sqcups: sqcups,
	Sqrt: Sqrt,
	sqsub: sqsub,
	sqsube: sqsube,
	sqsubset: sqsubset,
	sqsubseteq: sqsubseteq,
	sqsup: sqsup,
	sqsupe: sqsupe,
	sqsupset: sqsupset,
	sqsupseteq: sqsupseteq,
	square: square,
	Square: Square,
	SquareIntersection: SquareIntersection,
	SquareSubset: SquareSubset,
	SquareSubsetEqual: SquareSubsetEqual,
	SquareSuperset: SquareSuperset,
	SquareSupersetEqual: SquareSupersetEqual,
	SquareUnion: SquareUnion,
	squarf: squarf,
	squ: squ,
	squf: squf,
	srarr: srarr,
	Sscr: Sscr,
	sscr: sscr,
	ssetmn: ssetmn,
	ssmile: ssmile,
	sstarf: sstarf,
	Star: Star,
	star: star,
	starf: starf,
	straightepsilon: straightepsilon,
	straightphi: straightphi,
	strns: strns,
	sub: sub,
	Sub: Sub,
	subdot: subdot,
	subE: subE,
	sube: sube,
	subedot: subedot,
	submult: submult,
	subnE: subnE,
	subne: subne,
	subplus: subplus,
	subrarr: subrarr,
	subset: subset,
	Subset: Subset,
	subseteq: subseteq,
	subseteqq: subseteqq,
	SubsetEqual: SubsetEqual,
	subsetneq: subsetneq,
	subsetneqq: subsetneqq,
	subsim: subsim,
	subsub: subsub,
	subsup: subsup,
	succapprox: succapprox,
	succ: succ,
	succcurlyeq: succcurlyeq,
	Succeeds: Succeeds,
	SucceedsEqual: SucceedsEqual,
	SucceedsSlantEqual: SucceedsSlantEqual,
	SucceedsTilde: SucceedsTilde,
	succeq: succeq,
	succnapprox: succnapprox,
	succneqq: succneqq,
	succnsim: succnsim,
	succsim: succsim,
	SuchThat: SuchThat,
	sum: sum,
	Sum: Sum,
	sung: sung,
	sup1: sup1$1,
	sup2: sup2$1,
	sup3: sup3$1,
	sup: sup,
	Sup: Sup,
	supdot: supdot,
	supdsub: supdsub,
	supE: supE,
	supe: supe,
	supedot: supedot,
	Superset: Superset,
	SupersetEqual: SupersetEqual,
	suphsol: suphsol,
	suphsub: suphsub,
	suplarr: suplarr,
	supmult: supmult,
	supnE: supnE,
	supne: supne,
	supplus: supplus,
	supset: supset,
	Supset: Supset,
	supseteq: supseteq,
	supseteqq: supseteqq,
	supsetneq: supsetneq,
	supsetneqq: supsetneqq,
	supsim: supsim,
	supsub: supsub,
	supsup: supsup,
	swarhk: swarhk,
	swarr: swarr,
	swArr: swArr,
	swarrow: swarrow,
	swnwar: swnwar,
	szlig: szlig$1,
	Tab: Tab,
	target: target,
	Tau: Tau,
	tau: tau,
	tbrk: tbrk,
	Tcaron: Tcaron,
	tcaron: tcaron,
	Tcedil: Tcedil,
	tcedil: tcedil,
	Tcy: Tcy,
	tcy: tcy,
	tdot: tdot,
	telrec: telrec,
	Tfr: Tfr,
	tfr: tfr,
	there4: there4,
	therefore: therefore,
	Therefore: Therefore,
	Theta: Theta,
	theta: theta,
	thetasym: thetasym,
	thetav: thetav,
	thickapprox: thickapprox,
	thicksim: thicksim,
	ThickSpace: ThickSpace,
	ThinSpace: ThinSpace,
	thinsp: thinsp,
	thkap: thkap,
	thksim: thksim,
	THORN: THORN$1,
	thorn: thorn$1,
	tilde: tilde$1,
	Tilde: Tilde,
	TildeEqual: TildeEqual,
	TildeFullEqual: TildeFullEqual,
	TildeTilde: TildeTilde,
	timesbar: timesbar,
	timesb: timesb,
	times: times$1,
	timesd: timesd,
	tint: tint,
	toea: toea,
	topbot: topbot,
	topcir: topcir,
	top: top,
	Topf: Topf,
	topf: topf,
	topfork: topfork,
	tosa: tosa,
	tprime: tprime,
	trade: trade,
	TRADE: TRADE,
	triangle: triangle,
	triangledown: triangledown,
	triangleleft: triangleleft,
	trianglelefteq: trianglelefteq,
	triangleq: triangleq,
	triangleright: triangleright,
	trianglerighteq: trianglerighteq,
	tridot: tridot,
	trie: trie,
	triminus: triminus,
	TripleDot: TripleDot,
	triplus: triplus,
	trisb: trisb,
	tritime: tritime,
	trpezium: trpezium,
	Tscr: Tscr,
	tscr: tscr,
	TScy: TScy,
	tscy: tscy,
	TSHcy: TSHcy,
	tshcy: tshcy,
	Tstrok: Tstrok,
	tstrok: tstrok,
	twixt: twixt,
	twoheadleftarrow: twoheadleftarrow,
	twoheadrightarrow: twoheadrightarrow,
	Uacute: Uacute$1,
	uacute: uacute$1,
	uarr: uarr,
	Uarr: Uarr,
	uArr: uArr,
	Uarrocir: Uarrocir,
	Ubrcy: Ubrcy,
	ubrcy: ubrcy,
	Ubreve: Ubreve,
	ubreve: ubreve,
	Ucirc: Ucirc$1,
	ucirc: ucirc$1,
	Ucy: Ucy,
	ucy: ucy,
	udarr: udarr,
	Udblac: Udblac,
	udblac: udblac,
	udhar: udhar,
	ufisht: ufisht,
	Ufr: Ufr,
	ufr: ufr,
	Ugrave: Ugrave$1,
	ugrave: ugrave$1,
	uHar: uHar,
	uharl: uharl,
	uharr: uharr,
	uhblk: uhblk,
	ulcorn: ulcorn,
	ulcorner: ulcorner,
	ulcrop: ulcrop,
	ultri: ultri,
	Umacr: Umacr,
	umacr: umacr,
	uml: uml$1,
	UnderBar: UnderBar,
	UnderBrace: UnderBrace,
	UnderBracket: UnderBracket,
	UnderParenthesis: UnderParenthesis,
	Union: Union,
	UnionPlus: UnionPlus,
	Uogon: Uogon,
	uogon: uogon,
	Uopf: Uopf,
	uopf: uopf,
	UpArrowBar: UpArrowBar,
	uparrow: uparrow,
	UpArrow: UpArrow,
	Uparrow: Uparrow,
	UpArrowDownArrow: UpArrowDownArrow,
	updownarrow: updownarrow,
	UpDownArrow: UpDownArrow,
	Updownarrow: Updownarrow,
	UpEquilibrium: UpEquilibrium,
	upharpoonleft: upharpoonleft,
	upharpoonright: upharpoonright,
	uplus: uplus,
	UpperLeftArrow: UpperLeftArrow,
	UpperRightArrow: UpperRightArrow,
	upsi: upsi,
	Upsi: Upsi,
	upsih: upsih,
	Upsilon: Upsilon,
	upsilon: upsilon,
	UpTeeArrow: UpTeeArrow,
	UpTee: UpTee,
	upuparrows: upuparrows,
	urcorn: urcorn,
	urcorner: urcorner,
	urcrop: urcrop,
	Uring: Uring,
	uring: uring,
	urtri: urtri,
	Uscr: Uscr,
	uscr: uscr,
	utdot: utdot,
	Utilde: Utilde,
	utilde: utilde,
	utri: utri,
	utrif: utrif,
	uuarr: uuarr,
	Uuml: Uuml$1,
	uuml: uuml$1,
	uwangle: uwangle,
	vangrt: vangrt,
	varepsilon: varepsilon,
	varkappa: varkappa,
	varnothing: varnothing,
	varphi: varphi,
	varpi: varpi,
	varpropto: varpropto,
	varr: varr,
	vArr: vArr,
	varrho: varrho,
	varsigma: varsigma,
	varsubsetneq: varsubsetneq,
	varsubsetneqq: varsubsetneqq,
	varsupsetneq: varsupsetneq,
	varsupsetneqq: varsupsetneqq,
	vartheta: vartheta,
	vartriangleleft: vartriangleleft,
	vartriangleright: vartriangleright,
	vBar: vBar,
	Vbar: Vbar,
	vBarv: vBarv,
	Vcy: Vcy,
	vcy: vcy,
	vdash: vdash,
	vDash: vDash,
	Vdash: Vdash,
	VDash: VDash,
	Vdashl: Vdashl,
	veebar: veebar,
	vee: vee,
	Vee: Vee,
	veeeq: veeeq,
	vellip: vellip,
	verbar: verbar,
	Verbar: Verbar,
	vert: vert,
	Vert: Vert,
	VerticalBar: VerticalBar,
	VerticalLine: VerticalLine,
	VerticalSeparator: VerticalSeparator,
	VerticalTilde: VerticalTilde,
	VeryThinSpace: VeryThinSpace,
	Vfr: Vfr,
	vfr: vfr,
	vltri: vltri,
	vnsub: vnsub,
	vnsup: vnsup,
	Vopf: Vopf,
	vopf: vopf,
	vprop: vprop,
	vrtri: vrtri,
	Vscr: Vscr,
	vscr: vscr,
	vsubnE: vsubnE,
	vsubne: vsubne,
	vsupnE: vsupnE,
	vsupne: vsupne,
	Vvdash: Vvdash,
	vzigzag: vzigzag,
	Wcirc: Wcirc,
	wcirc: wcirc,
	wedbar: wedbar,
	wedge: wedge,
	Wedge: Wedge,
	wedgeq: wedgeq,
	weierp: weierp,
	Wfr: Wfr,
	wfr: wfr,
	Wopf: Wopf,
	wopf: wopf,
	wp: wp,
	wr: wr,
	wreath: wreath,
	Wscr: Wscr,
	wscr: wscr,
	xcap: xcap,
	xcirc: xcirc,
	xcup: xcup,
	xdtri: xdtri,
	Xfr: Xfr,
	xfr: xfr,
	xharr: xharr,
	xhArr: xhArr,
	Xi: Xi,
	xi: xi,
	xlarr: xlarr,
	xlArr: xlArr,
	xmap: xmap,
	xnis: xnis,
	xodot: xodot,
	Xopf: Xopf,
	xopf: xopf,
	xoplus: xoplus,
	xotime: xotime,
	xrarr: xrarr,
	xrArr: xrArr,
	Xscr: Xscr,
	xscr: xscr,
	xsqcup: xsqcup,
	xuplus: xuplus,
	xutri: xutri,
	xvee: xvee,
	xwedge: xwedge,
	Yacute: Yacute$1,
	yacute: yacute$1,
	YAcy: YAcy,
	yacy: yacy,
	Ycirc: Ycirc,
	ycirc: ycirc,
	Ycy: Ycy,
	ycy: ycy,
	yen: yen$1,
	Yfr: Yfr,
	yfr: yfr,
	YIcy: YIcy,
	yicy: yicy,
	Yopf: Yopf,
	yopf: yopf,
	Yscr: Yscr,
	yscr: yscr,
	YUcy: YUcy,
	yucy: yucy,
	yuml: yuml$1,
	Yuml: Yuml,
	Zacute: Zacute,
	zacute: zacute,
	Zcaron: Zcaron,
	zcaron: zcaron,
	Zcy: Zcy,
	zcy: zcy,
	Zdot: Zdot,
	zdot: zdot,
	zeetrf: zeetrf,
	ZeroWidthSpace: ZeroWidthSpace,
	Zeta: Zeta,
	zeta: zeta,
	zfr: zfr,
	Zfr: Zfr,
	ZHcy: ZHcy,
	zhcy: zhcy,
	zigrarr: zigrarr,
	zopf: zopf,
	Zopf: Zopf,
	Zscr: Zscr,
	zscr: zscr,
	zwj: zwj,
	zwnj: zwnj
};

var Aacute = "Á";
var aacute = "á";
var Acirc = "Â";
var acirc = "â";
var acute = "´";
var AElig = "Æ";
var aelig = "æ";
var Agrave = "À";
var agrave = "à";
var amp$1 = "&";
var AMP = "&";
var Aring = "Å";
var aring = "å";
var Atilde = "Ã";
var atilde = "ã";
var Auml = "Ä";
var auml = "ä";
var brvbar = "¦";
var Ccedil = "Ç";
var ccedil = "ç";
var cedil = "¸";
var cent = "¢";
var copy = "©";
var COPY = "©";
var curren = "¤";
var deg = "°";
var divide = "÷";
var Eacute = "É";
var eacute = "é";
var Ecirc = "Ê";
var ecirc = "ê";
var Egrave = "È";
var egrave = "è";
var ETH = "Ð";
var eth = "ð";
var Euml = "Ë";
var euml = "ë";
var frac12 = "½";
var frac14 = "¼";
var frac34 = "¾";
var gt$2 = ">";
var GT = ">";
var Iacute = "Í";
var iacute = "í";
var Icirc = "Î";
var icirc = "î";
var iexcl = "¡";
var Igrave = "Ì";
var igrave = "ì";
var iquest = "¿";
var Iuml = "Ï";
var iuml = "ï";
var laquo = "«";
var lt$1 = "<";
var LT = "<";
var macr = "¯";
var micro = "µ";
var middot = "·";
var nbsp = " ";
var not = "¬";
var Ntilde = "Ñ";
var ntilde = "ñ";
var Oacute = "Ó";
var oacute = "ó";
var Ocirc = "Ô";
var ocirc = "ô";
var Ograve = "Ò";
var ograve = "ò";
var ordf = "ª";
var ordm = "º";
var Oslash = "Ø";
var oslash = "ø";
var Otilde = "Õ";
var otilde = "õ";
var Ouml = "Ö";
var ouml = "ö";
var para = "¶";
var plusmn = "±";
var pound = "£";
var quot$1 = "\"";
var QUOT = "\"";
var raquo = "»";
var reg = "®";
var REG = "®";
var sect = "§";
var shy = "­";
var sup1 = "¹";
var sup2 = "²";
var sup3 = "³";
var szlig = "ß";
var THORN = "Þ";
var thorn = "þ";
var times = "×";
var Uacute = "Ú";
var uacute = "ú";
var Ucirc = "Û";
var ucirc = "û";
var Ugrave = "Ù";
var ugrave = "ù";
var uml = "¨";
var Uuml = "Ü";
var uuml = "ü";
var Yacute = "Ý";
var yacute = "ý";
var yen = "¥";
var yuml = "ÿ";
var require$$1 = {
	Aacute: Aacute,
	aacute: aacute,
	Acirc: Acirc,
	acirc: acirc,
	acute: acute,
	AElig: AElig,
	aelig: aelig,
	Agrave: Agrave,
	agrave: agrave,
	amp: amp$1,
	AMP: AMP,
	Aring: Aring,
	aring: aring,
	Atilde: Atilde,
	atilde: atilde,
	Auml: Auml,
	auml: auml,
	brvbar: brvbar,
	Ccedil: Ccedil,
	ccedil: ccedil,
	cedil: cedil,
	cent: cent,
	copy: copy,
	COPY: COPY,
	curren: curren,
	deg: deg,
	divide: divide,
	Eacute: Eacute,
	eacute: eacute,
	Ecirc: Ecirc,
	ecirc: ecirc,
	Egrave: Egrave,
	egrave: egrave,
	ETH: ETH,
	eth: eth,
	Euml: Euml,
	euml: euml,
	frac12: frac12,
	frac14: frac14,
	frac34: frac34,
	gt: gt$2,
	GT: GT,
	Iacute: Iacute,
	iacute: iacute,
	Icirc: Icirc,
	icirc: icirc,
	iexcl: iexcl,
	Igrave: Igrave,
	igrave: igrave,
	iquest: iquest,
	Iuml: Iuml,
	iuml: iuml,
	laquo: laquo,
	lt: lt$1,
	LT: LT,
	macr: macr,
	micro: micro,
	middot: middot,
	nbsp: nbsp,
	not: not,
	Ntilde: Ntilde,
	ntilde: ntilde,
	Oacute: Oacute,
	oacute: oacute,
	Ocirc: Ocirc,
	ocirc: ocirc,
	Ograve: Ograve,
	ograve: ograve,
	ordf: ordf,
	ordm: ordm,
	Oslash: Oslash,
	oslash: oslash,
	Otilde: Otilde,
	otilde: otilde,
	Ouml: Ouml,
	ouml: ouml,
	para: para,
	plusmn: plusmn,
	pound: pound,
	quot: quot$1,
	QUOT: QUOT,
	raquo: raquo,
	reg: reg,
	REG: REG,
	sect: sect,
	shy: shy,
	sup1: sup1,
	sup2: sup2,
	sup3: sup3,
	szlig: szlig,
	THORN: THORN,
	thorn: thorn,
	times: times,
	Uacute: Uacute,
	uacute: uacute,
	Ucirc: Ucirc,
	ucirc: ucirc,
	Ugrave: Ugrave,
	ugrave: ugrave,
	uml: uml,
	Uuml: Uuml,
	uuml: uuml,
	Yacute: Yacute,
	yacute: yacute,
	yen: yen,
	yuml: yuml
};

var amp = "&";
var apos = "'";
var gt$1 = ">";
var lt = "<";
var quot = "\"";
var require$$0$1 = {
	amp: amp,
	apos: apos,
	gt: gt$1,
	lt: lt,
	quot: quot
};

var decode_codepoint = {};

var require$$0 = {
	"0": 65533,
	"128": 8364,
	"130": 8218,
	"131": 402,
	"132": 8222,
	"133": 8230,
	"134": 8224,
	"135": 8225,
	"136": 710,
	"137": 8240,
	"138": 352,
	"139": 8249,
	"140": 338,
	"142": 381,
	"145": 8216,
	"146": 8217,
	"147": 8220,
	"148": 8221,
	"149": 8226,
	"150": 8211,
	"151": 8212,
	"152": 732,
	"153": 8482,
	"154": 353,
	"155": 8250,
	"156": 339,
	"158": 382,
	"159": 376
};

var __importDefault$4 = (commonjsGlobal && commonjsGlobal.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(decode_codepoint, "__esModule", { value: true });
var decode_json_1 = __importDefault$4(require$$0);
// Adapted from https://github.com/mathiasbynens/he/blob/master/src/he.js#L94-L119
var fromCodePoint = 
// eslint-disable-next-line @typescript-eslint/no-unnecessary-condition
String.fromCodePoint ||
    function (codePoint) {
        var output = "";
        if (codePoint > 0xffff) {
            codePoint -= 0x10000;
            output += String.fromCharCode(((codePoint >>> 10) & 0x3ff) | 0xd800);
            codePoint = 0xdc00 | (codePoint & 0x3ff);
        }
        output += String.fromCharCode(codePoint);
        return output;
    };
function decodeCodePoint(codePoint) {
    if ((codePoint >= 0xd800 && codePoint <= 0xdfff) || codePoint > 0x10ffff) {
        return "\uFFFD";
    }
    if (codePoint in decode_json_1.default) {
        codePoint = decode_json_1.default[codePoint];
    }
    return fromCodePoint(codePoint);
}
decode_codepoint.default = decodeCodePoint;

var __importDefault$3 = (commonjsGlobal && commonjsGlobal.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(decode, "__esModule", { value: true });
decode.decodeHTML = decode.decodeHTMLStrict = decode.decodeXML = void 0;
var entities_json_1$1 = __importDefault$3(require$$1$1);
var legacy_json_1 = __importDefault$3(require$$1);
var xml_json_1$1 = __importDefault$3(require$$0$1);
var decode_codepoint_1 = __importDefault$3(decode_codepoint);
var strictEntityRe = /&(?:[a-zA-Z0-9]+|#[xX][\da-fA-F]+|#\d+);/g;
decode.decodeXML = getStrictDecoder(xml_json_1$1.default);
decode.decodeHTMLStrict = getStrictDecoder(entities_json_1$1.default);
function getStrictDecoder(map) {
    var replace = getReplacer(map);
    return function (str) { return String(str).replace(strictEntityRe, replace); };
}
var sorter = function (a, b) { return (a < b ? 1 : -1); };
decode.decodeHTML = (function () {
    var legacy = Object.keys(legacy_json_1.default).sort(sorter);
    var keys = Object.keys(entities_json_1$1.default).sort(sorter);
    for (var i = 0, j = 0; i < keys.length; i++) {
        if (legacy[j] === keys[i]) {
            keys[i] += ";?";
            j++;
        }
        else {
            keys[i] += ";";
        }
    }
    var re = new RegExp("&(?:" + keys.join("|") + "|#[xX][\\da-fA-F]+;?|#\\d+;?)", "g");
    var replace = getReplacer(entities_json_1$1.default);
    function replacer(str) {
        if (str.substr(-1) !== ";")
            str += ";";
        return replace(str);
    }
    // TODO consider creating a merged map
    return function (str) { return String(str).replace(re, replacer); };
})();
function getReplacer(map) {
    return function replace(str) {
        if (str.charAt(1) === "#") {
            var secondChar = str.charAt(2);
            if (secondChar === "X" || secondChar === "x") {
                return decode_codepoint_1.default(parseInt(str.substr(3), 16));
            }
            return decode_codepoint_1.default(parseInt(str.substr(2), 10));
        }
        // eslint-disable-next-line @typescript-eslint/prefer-nullish-coalescing
        return map[str.slice(1, -1)] || str;
    };
}

var encode = {};

var __importDefault$2 = (commonjsGlobal && commonjsGlobal.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(encode, "__esModule", { value: true });
encode.escapeUTF8 = encode.escape = encode.encodeNonAsciiHTML = encode.encodeHTML = encode.encodeXML = void 0;
var xml_json_1 = __importDefault$2(require$$0$1);
var inverseXML = getInverseObj(xml_json_1.default);
var xmlReplacer = getInverseReplacer(inverseXML);
/**
 * Encodes all non-ASCII characters, as well as characters not valid in XML
 * documents using XML entities.
 *
 * If a character has no equivalent entity, a
 * numeric hexadecimal reference (eg. `&#xfc;`) will be used.
 */
encode.encodeXML = getASCIIEncoder(inverseXML);
var entities_json_1 = __importDefault$2(require$$1$1);
var inverseHTML = getInverseObj(entities_json_1.default);
var htmlReplacer = getInverseReplacer(inverseHTML);
/**
 * Encodes all entities and non-ASCII characters in the input.
 *
 * This includes characters that are valid ASCII characters in HTML documents.
 * For example `#` will be encoded as `&num;`. To get a more compact output,
 * consider using the `encodeNonAsciiHTML` function.
 *
 * If a character has no equivalent entity, a
 * numeric hexadecimal reference (eg. `&#xfc;`) will be used.
 */
encode.encodeHTML = getInverse(inverseHTML, htmlReplacer);
/**
 * Encodes all non-ASCII characters, as well as characters not valid in HTML
 * documents using HTML entities.
 *
 * If a character has no equivalent entity, a
 * numeric hexadecimal reference (eg. `&#xfc;`) will be used.
 */
encode.encodeNonAsciiHTML = getASCIIEncoder(inverseHTML);
function getInverseObj(obj) {
    return Object.keys(obj)
        .sort()
        .reduce(function (inverse, name) {
        inverse[obj[name]] = "&" + name + ";";
        return inverse;
    }, {});
}
function getInverseReplacer(inverse) {
    var single = [];
    var multiple = [];
    for (var _i = 0, _a = Object.keys(inverse); _i < _a.length; _i++) {
        var k = _a[_i];
        if (k.length === 1) {
            // Add value to single array
            single.push("\\" + k);
        }
        else {
            // Add value to multiple array
            multiple.push(k);
        }
    }
    // Add ranges to single characters.
    single.sort();
    for (var start = 0; start < single.length - 1; start++) {
        // Find the end of a run of characters
        var end = start;
        while (end < single.length - 1 &&
            single[end].charCodeAt(1) + 1 === single[end + 1].charCodeAt(1)) {
            end += 1;
        }
        var count = 1 + end - start;
        // We want to replace at least three characters
        if (count < 3)
            continue;
        single.splice(start, count, single[start] + "-" + single[end]);
    }
    multiple.unshift("[" + single.join("") + "]");
    return new RegExp(multiple.join("|"), "g");
}
// /[^\0-\x7F]/gu
var reNonASCII = /(?:[\x80-\uD7FF\uE000-\uFFFF]|[\uD800-\uDBFF][\uDC00-\uDFFF]|[\uD800-\uDBFF](?![\uDC00-\uDFFF])|(?:[^\uD800-\uDBFF]|^)[\uDC00-\uDFFF])/g;
var getCodePoint = 
// eslint-disable-next-line @typescript-eslint/no-unnecessary-condition
String.prototype.codePointAt != null
    ? // eslint-disable-next-line @typescript-eslint/no-non-null-assertion
        function (str) { return str.codePointAt(0); }
    : // http://mathiasbynens.be/notes/javascript-encoding#surrogate-formulae
        function (c) {
            return (c.charCodeAt(0) - 0xd800) * 0x400 +
                c.charCodeAt(1) -
                0xdc00 +
                0x10000;
        };
function singleCharReplacer(c) {
    return "&#x" + (c.length > 1 ? getCodePoint(c) : c.charCodeAt(0))
        .toString(16)
        .toUpperCase() + ";";
}
function getInverse(inverse, re) {
    return function (data) {
        return data
            .replace(re, function (name) { return inverse[name]; })
            .replace(reNonASCII, singleCharReplacer);
    };
}
var reEscapeChars = new RegExp(xmlReplacer.source + "|" + reNonASCII.source, "g");
/**
 * Encodes all non-ASCII characters, as well as characters not valid in XML
 * documents using numeric hexadecimal reference (eg. `&#xfc;`).
 *
 * Have a look at `escapeUTF8` if you want a more concise output at the expense
 * of reduced transportability.
 *
 * @param data String to escape.
 */
function escape(data) {
    return data.replace(reEscapeChars, singleCharReplacer);
}
encode.escape = escape;
/**
 * Encodes all characters not valid in XML documents using numeric hexadecimal
 * reference (eg. `&#xfc;`).
 *
 * Note that the output will be character-set dependent.
 *
 * @param data String to escape.
 */
function escapeUTF8(data) {
    return data.replace(xmlReplacer, singleCharReplacer);
}
encode.escapeUTF8 = escapeUTF8;
function getASCIIEncoder(obj) {
    return function (data) {
        return data.replace(reEscapeChars, function (c) { return obj[c] || singleCharReplacer(c); });
    };
}

(function (exports) {
	Object.defineProperty(exports, "__esModule", { value: true });
	exports.decodeXMLStrict = exports.decodeHTML5Strict = exports.decodeHTML4Strict = exports.decodeHTML5 = exports.decodeHTML4 = exports.decodeHTMLStrict = exports.decodeHTML = exports.decodeXML = exports.encodeHTML5 = exports.encodeHTML4 = exports.escapeUTF8 = exports.escape = exports.encodeNonAsciiHTML = exports.encodeHTML = exports.encodeXML = exports.encode = exports.decodeStrict = exports.decode = void 0;
	var decode_1 = decode;
	var encode_1 = encode;
	/**
	 * Decodes a string with entities.
	 *
	 * @param data String to decode.
	 * @param level Optional level to decode at. 0 = XML, 1 = HTML. Default is 0.
	 * @deprecated Use `decodeXML` or `decodeHTML` directly.
	 */
	function decode$1(data, level) {
	    return (!level || level <= 0 ? decode_1.decodeXML : decode_1.decodeHTML)(data);
	}
	exports.decode = decode$1;
	/**
	 * Decodes a string with entities. Does not allow missing trailing semicolons for entities.
	 *
	 * @param data String to decode.
	 * @param level Optional level to decode at. 0 = XML, 1 = HTML. Default is 0.
	 * @deprecated Use `decodeHTMLStrict` or `decodeXML` directly.
	 */
	function decodeStrict(data, level) {
	    return (!level || level <= 0 ? decode_1.decodeXML : decode_1.decodeHTMLStrict)(data);
	}
	exports.decodeStrict = decodeStrict;
	/**
	 * Encodes a string with entities.
	 *
	 * @param data String to encode.
	 * @param level Optional level to encode at. 0 = XML, 1 = HTML. Default is 0.
	 * @deprecated Use `encodeHTML`, `encodeXML` or `encodeNonAsciiHTML` directly.
	 */
	function encode$1(data, level) {
	    return (!level || level <= 0 ? encode_1.encodeXML : encode_1.encodeHTML)(data);
	}
	exports.encode = encode$1;
	var encode_2 = encode;
	Object.defineProperty(exports, "encodeXML", { enumerable: true, get: function () { return encode_2.encodeXML; } });
	Object.defineProperty(exports, "encodeHTML", { enumerable: true, get: function () { return encode_2.encodeHTML; } });
	Object.defineProperty(exports, "encodeNonAsciiHTML", { enumerable: true, get: function () { return encode_2.encodeNonAsciiHTML; } });
	Object.defineProperty(exports, "escape", { enumerable: true, get: function () { return encode_2.escape; } });
	Object.defineProperty(exports, "escapeUTF8", { enumerable: true, get: function () { return encode_2.escapeUTF8; } });
	// Legacy aliases (deprecated)
	Object.defineProperty(exports, "encodeHTML4", { enumerable: true, get: function () { return encode_2.encodeHTML; } });
	Object.defineProperty(exports, "encodeHTML5", { enumerable: true, get: function () { return encode_2.encodeHTML; } });
	var decode_2 = decode;
	Object.defineProperty(exports, "decodeXML", { enumerable: true, get: function () { return decode_2.decodeXML; } });
	Object.defineProperty(exports, "decodeHTML", { enumerable: true, get: function () { return decode_2.decodeHTML; } });
	Object.defineProperty(exports, "decodeHTMLStrict", { enumerable: true, get: function () { return decode_2.decodeHTMLStrict; } });
	// Legacy aliases (deprecated)
	Object.defineProperty(exports, "decodeHTML4", { enumerable: true, get: function () { return decode_2.decodeHTML; } });
	Object.defineProperty(exports, "decodeHTML5", { enumerable: true, get: function () { return decode_2.decodeHTML; } });
	Object.defineProperty(exports, "decodeHTML4Strict", { enumerable: true, get: function () { return decode_2.decodeHTMLStrict; } });
	Object.defineProperty(exports, "decodeHTML5Strict", { enumerable: true, get: function () { return decode_2.decodeHTMLStrict; } });
	Object.defineProperty(exports, "decodeXMLStrict", { enumerable: true, get: function () { return decode_2.decodeXML; } }); 
} (lib));

var foreignNames = {};

Object.defineProperty(foreignNames, "__esModule", { value: true });
foreignNames.attributeNames = foreignNames.elementNames = void 0;
foreignNames.elementNames = new Map([
    ["altglyph", "altGlyph"],
    ["altglyphdef", "altGlyphDef"],
    ["altglyphitem", "altGlyphItem"],
    ["animatecolor", "animateColor"],
    ["animatemotion", "animateMotion"],
    ["animatetransform", "animateTransform"],
    ["clippath", "clipPath"],
    ["feblend", "feBlend"],
    ["fecolormatrix", "feColorMatrix"],
    ["fecomponenttransfer", "feComponentTransfer"],
    ["fecomposite", "feComposite"],
    ["feconvolvematrix", "feConvolveMatrix"],
    ["fediffuselighting", "feDiffuseLighting"],
    ["fedisplacementmap", "feDisplacementMap"],
    ["fedistantlight", "feDistantLight"],
    ["fedropshadow", "feDropShadow"],
    ["feflood", "feFlood"],
    ["fefunca", "feFuncA"],
    ["fefuncb", "feFuncB"],
    ["fefuncg", "feFuncG"],
    ["fefuncr", "feFuncR"],
    ["fegaussianblur", "feGaussianBlur"],
    ["feimage", "feImage"],
    ["femerge", "feMerge"],
    ["femergenode", "feMergeNode"],
    ["femorphology", "feMorphology"],
    ["feoffset", "feOffset"],
    ["fepointlight", "fePointLight"],
    ["fespecularlighting", "feSpecularLighting"],
    ["fespotlight", "feSpotLight"],
    ["fetile", "feTile"],
    ["feturbulence", "feTurbulence"],
    ["foreignobject", "foreignObject"],
    ["glyphref", "glyphRef"],
    ["lineargradient", "linearGradient"],
    ["radialgradient", "radialGradient"],
    ["textpath", "textPath"],
]);
foreignNames.attributeNames = new Map([
    ["definitionurl", "definitionURL"],
    ["attributename", "attributeName"],
    ["attributetype", "attributeType"],
    ["basefrequency", "baseFrequency"],
    ["baseprofile", "baseProfile"],
    ["calcmode", "calcMode"],
    ["clippathunits", "clipPathUnits"],
    ["diffuseconstant", "diffuseConstant"],
    ["edgemode", "edgeMode"],
    ["filterunits", "filterUnits"],
    ["glyphref", "glyphRef"],
    ["gradienttransform", "gradientTransform"],
    ["gradientunits", "gradientUnits"],
    ["kernelmatrix", "kernelMatrix"],
    ["kernelunitlength", "kernelUnitLength"],
    ["keypoints", "keyPoints"],
    ["keysplines", "keySplines"],
    ["keytimes", "keyTimes"],
    ["lengthadjust", "lengthAdjust"],
    ["limitingconeangle", "limitingConeAngle"],
    ["markerheight", "markerHeight"],
    ["markerunits", "markerUnits"],
    ["markerwidth", "markerWidth"],
    ["maskcontentunits", "maskContentUnits"],
    ["maskunits", "maskUnits"],
    ["numoctaves", "numOctaves"],
    ["pathlength", "pathLength"],
    ["patterncontentunits", "patternContentUnits"],
    ["patterntransform", "patternTransform"],
    ["patternunits", "patternUnits"],
    ["pointsatx", "pointsAtX"],
    ["pointsaty", "pointsAtY"],
    ["pointsatz", "pointsAtZ"],
    ["preservealpha", "preserveAlpha"],
    ["preserveaspectratio", "preserveAspectRatio"],
    ["primitiveunits", "primitiveUnits"],
    ["refx", "refX"],
    ["refy", "refY"],
    ["repeatcount", "repeatCount"],
    ["repeatdur", "repeatDur"],
    ["requiredextensions", "requiredExtensions"],
    ["requiredfeatures", "requiredFeatures"],
    ["specularconstant", "specularConstant"],
    ["specularexponent", "specularExponent"],
    ["spreadmethod", "spreadMethod"],
    ["startoffset", "startOffset"],
    ["stddeviation", "stdDeviation"],
    ["stitchtiles", "stitchTiles"],
    ["surfacescale", "surfaceScale"],
    ["systemlanguage", "systemLanguage"],
    ["tablevalues", "tableValues"],
    ["targetx", "targetX"],
    ["targety", "targetY"],
    ["textlength", "textLength"],
    ["viewbox", "viewBox"],
    ["viewtarget", "viewTarget"],
    ["xchannelselector", "xChannelSelector"],
    ["ychannelselector", "yChannelSelector"],
    ["zoomandpan", "zoomAndPan"],
]);

var __assign = (commonjsGlobal && commonjsGlobal.__assign) || function () {
    __assign = Object.assign || function(t) {
        for (var s, i = 1, n = arguments.length; i < n; i++) {
            s = arguments[i];
            for (var p in s) if (Object.prototype.hasOwnProperty.call(s, p))
                t[p] = s[p];
        }
        return t;
    };
    return __assign.apply(this, arguments);
};
var __createBinding$1 = (commonjsGlobal && commonjsGlobal.__createBinding) || (Object.create ? (function(o, m, k, k2) {
    if (k2 === undefined) k2 = k;
    Object.defineProperty(o, k2, { enumerable: true, get: function() { return m[k]; } });
}) : (function(o, m, k, k2) {
    if (k2 === undefined) k2 = k;
    o[k2] = m[k];
}));
var __setModuleDefault$1 = (commonjsGlobal && commonjsGlobal.__setModuleDefault) || (Object.create ? (function(o, v) {
    Object.defineProperty(o, "default", { enumerable: true, value: v });
}) : function(o, v) {
    o["default"] = v;
});
var __importStar$1 = (commonjsGlobal && commonjsGlobal.__importStar) || function (mod) {
    if (mod && mod.__esModule) return mod;
    var result = {};
    if (mod != null) for (var k in mod) if (k !== "default" && Object.prototype.hasOwnProperty.call(mod, k)) __createBinding$1(result, mod, k);
    __setModuleDefault$1(result, mod);
    return result;
};
Object.defineProperty(lib$1, "__esModule", { value: true });
/*
 * Module dependencies
 */
var ElementType = __importStar$1(lib$4);
var entities_1 = lib;
/**
 * Mixed-case SVG and MathML tags & attributes
 * recognized by the HTML parser.
 *
 * @see https://html.spec.whatwg.org/multipage/parsing.html#parsing-main-inforeign
 */
var foreignNames_1 = foreignNames;
var unencodedElements = new Set([
    "style",
    "script",
    "xmp",
    "iframe",
    "noembed",
    "noframes",
    "plaintext",
    "noscript",
]);
/**
 * Format attributes
 */
function formatAttributes(attributes, opts) {
    if (!attributes)
        return;
    return Object.keys(attributes)
        .map(function (key) {
        var _a, _b;
        var value = (_a = attributes[key]) !== null && _a !== void 0 ? _a : "";
        if (opts.xmlMode === "foreign") {
            /* Fix up mixed-case attribute names */
            key = (_b = foreignNames_1.attributeNames.get(key)) !== null && _b !== void 0 ? _b : key;
        }
        if (!opts.emptyAttrs && !opts.xmlMode && value === "") {
            return key;
        }
        return key + "=\"" + (opts.decodeEntities !== false
            ? entities_1.encodeXML(value)
            : value.replace(/"/g, "&quot;")) + "\"";
    })
        .join(" ");
}
/**
 * Self-enclosing tags
 */
var singleTag = new Set([
    "area",
    "base",
    "basefont",
    "br",
    "col",
    "command",
    "embed",
    "frame",
    "hr",
    "img",
    "input",
    "isindex",
    "keygen",
    "link",
    "meta",
    "param",
    "source",
    "track",
    "wbr",
]);
/**
 * Renders a DOM node or an array of DOM nodes to a string.
 *
 * Can be thought of as the equivalent of the `outerHTML` of the passed node(s).
 *
 * @param node Node to be rendered.
 * @param options Changes serialization behavior
 */
function render(node, options) {
    if (options === void 0) { options = {}; }
    var nodes = "length" in node ? node : [node];
    var output = "";
    for (var i = 0; i < nodes.length; i++) {
        output += renderNode(nodes[i], options);
    }
    return output;
}
lib$1.default = render;
function renderNode(node, options) {
    switch (node.type) {
        case ElementType.Root:
            return render(node.children, options);
        case ElementType.Directive:
        case ElementType.Doctype:
            return renderDirective(node);
        case ElementType.Comment:
            return renderComment(node);
        case ElementType.CDATA:
            return renderCdata(node);
        case ElementType.Script:
        case ElementType.Style:
        case ElementType.Tag:
            return renderTag(node, options);
        case ElementType.Text:
            return renderText(node, options);
    }
}
var foreignModeIntegrationPoints = new Set([
    "mi",
    "mo",
    "mn",
    "ms",
    "mtext",
    "annotation-xml",
    "foreignObject",
    "desc",
    "title",
]);
var foreignElements = new Set(["svg", "math"]);
function renderTag(elem, opts) {
    var _a;
    // Handle SVG / MathML in HTML
    if (opts.xmlMode === "foreign") {
        /* Fix up mixed-case element names */
        elem.name = (_a = foreignNames_1.elementNames.get(elem.name)) !== null && _a !== void 0 ? _a : elem.name;
        /* Exit foreign mode at integration points */
        if (elem.parent &&
            foreignModeIntegrationPoints.has(elem.parent.name)) {
            opts = __assign(__assign({}, opts), { xmlMode: false });
        }
    }
    if (!opts.xmlMode && foreignElements.has(elem.name)) {
        opts = __assign(__assign({}, opts), { xmlMode: "foreign" });
    }
    var tag = "<" + elem.name;
    var attribs = formatAttributes(elem.attribs, opts);
    if (attribs) {
        tag += " " + attribs;
    }
    if (elem.children.length === 0 &&
        (opts.xmlMode
            ? // In XML mode or foreign mode, and user hasn't explicitly turned off self-closing tags
                opts.selfClosingTags !== false
            : // User explicitly asked for self-closing tags, even in HTML mode
                opts.selfClosingTags && singleTag.has(elem.name))) {
        if (!opts.xmlMode)
            tag += " ";
        tag += "/>";
    }
    else {
        tag += ">";
        if (elem.children.length > 0) {
            tag += render(elem.children, opts);
        }
        if (opts.xmlMode || !singleTag.has(elem.name)) {
            tag += "</" + elem.name + ">";
        }
    }
    return tag;
}
function renderDirective(elem) {
    return "<" + elem.data + ">";
}
function renderText(elem, opts) {
    var data = elem.data || "";
    // If entities weren't decoded, no need to encode them back
    if (opts.decodeEntities !== false &&
        !(!opts.xmlMode &&
            elem.parent &&
            unencodedElements.has(elem.parent.name))) {
        data = entities_1.encodeXML(data);
    }
    return data;
}
function renderCdata(elem) {
    return "<![CDATA[" + elem.children[0].data + "]]>";
}
function renderComment(elem) {
    return "<!--" + elem.data + "-->";
}

var __importDefault$1 = (commonjsGlobal && commonjsGlobal.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(stringify, "__esModule", { value: true });
stringify.innerText = stringify.textContent = stringify.getText = stringify.getInnerHTML = stringify.getOuterHTML = void 0;
var domhandler_1$5 = lib$5;
var dom_serializer_1 = __importDefault$1(lib$1);
var domelementtype_1 = lib$4;
/**
 * @param node Node to get the outer HTML of.
 * @param options Options for serialization.
 * @deprecated Use the `dom-serializer` module directly.
 * @returns `node`'s outer HTML.
 */
function getOuterHTML(node, options) {
    return (0, dom_serializer_1.default)(node, options);
}
stringify.getOuterHTML = getOuterHTML;
/**
 * @param node Node to get the inner HTML of.
 * @param options Options for serialization.
 * @deprecated Use the `dom-serializer` module directly.
 * @returns `node`'s inner HTML.
 */
function getInnerHTML(node, options) {
    return (0, domhandler_1$5.hasChildren)(node)
        ? node.children.map(function (node) { return getOuterHTML(node, options); }).join("")
        : "";
}
stringify.getInnerHTML = getInnerHTML;
/**
 * Get a node's inner text. Same as `textContent`, but inserts newlines for `<br>` tags.
 *
 * @deprecated Use `textContent` instead.
 * @param node Node to get the inner text of.
 * @returns `node`'s inner text.
 */
function getText$1(node) {
    if (Array.isArray(node))
        return node.map(getText$1).join("");
    if ((0, domhandler_1$5.isTag)(node))
        return node.name === "br" ? "\n" : getText$1(node.children);
    if ((0, domhandler_1$5.isCDATA)(node))
        return getText$1(node.children);
    if ((0, domhandler_1$5.isText)(node))
        return node.data;
    return "";
}
stringify.getText = getText$1;
/**
 * Get a node's text content.
 *
 * @param node Node to get the text content of.
 * @returns `node`'s text content.
 * @see {@link https://developer.mozilla.org/en-US/docs/Web/API/Node/textContent}
 */
function textContent(node) {
    if (Array.isArray(node))
        return node.map(textContent).join("");
    if ((0, domhandler_1$5.hasChildren)(node) && !(0, domhandler_1$5.isComment)(node)) {
        return textContent(node.children);
    }
    if ((0, domhandler_1$5.isText)(node))
        return node.data;
    return "";
}
stringify.textContent = textContent;
/**
 * Get a node's inner text.
 *
 * @param node Node to get the inner text of.
 * @returns `node`'s inner text.
 * @see {@link https://developer.mozilla.org/en-US/docs/Web/API/Node/innerText}
 */
function innerText(node) {
    if (Array.isArray(node))
        return node.map(innerText).join("");
    if ((0, domhandler_1$5.hasChildren)(node) && (node.type === domelementtype_1.ElementType.Tag || (0, domhandler_1$5.isCDATA)(node))) {
        return innerText(node.children);
    }
    if ((0, domhandler_1$5.isText)(node))
        return node.data;
    return "";
}
stringify.innerText = innerText;

var traversal = {};

Object.defineProperty(traversal, "__esModule", { value: true });
traversal.prevElementSibling = traversal.nextElementSibling = traversal.getName = traversal.hasAttrib = traversal.getAttributeValue = traversal.getSiblings = traversal.getParent = traversal.getChildren = void 0;
var domhandler_1$4 = lib$5;
var emptyArray = [];
/**
 * Get a node's children.
 *
 * @param elem Node to get the children of.
 * @returns `elem`'s children, or an empty array.
 */
function getChildren(elem) {
    var _a;
    return (_a = elem.children) !== null && _a !== void 0 ? _a : emptyArray;
}
traversal.getChildren = getChildren;
/**
 * Get a node's parent.
 *
 * @param elem Node to get the parent of.
 * @returns `elem`'s parent node.
 */
function getParent(elem) {
    return elem.parent || null;
}
traversal.getParent = getParent;
/**
 * Gets an elements siblings, including the element itself.
 *
 * Attempts to get the children through the element's parent first.
 * If we don't have a parent (the element is a root node),
 * we walk the element's `prev` & `next` to get all remaining nodes.
 *
 * @param elem Element to get the siblings of.
 * @returns `elem`'s siblings.
 */
function getSiblings(elem) {
    var _a, _b;
    var parent = getParent(elem);
    if (parent != null)
        return getChildren(parent);
    var siblings = [elem];
    var prev = elem.prev, next = elem.next;
    while (prev != null) {
        siblings.unshift(prev);
        (_a = prev, prev = _a.prev);
    }
    while (next != null) {
        siblings.push(next);
        (_b = next, next = _b.next);
    }
    return siblings;
}
traversal.getSiblings = getSiblings;
/**
 * Gets an attribute from an element.
 *
 * @param elem Element to check.
 * @param name Attribute name to retrieve.
 * @returns The element's attribute value, or `undefined`.
 */
function getAttributeValue(elem, name) {
    var _a;
    return (_a = elem.attribs) === null || _a === void 0 ? void 0 : _a[name];
}
traversal.getAttributeValue = getAttributeValue;
/**
 * Checks whether an element has an attribute.
 *
 * @param elem Element to check.
 * @param name Attribute name to look for.
 * @returns Returns whether `elem` has the attribute `name`.
 */
function hasAttrib(elem, name) {
    return (elem.attribs != null &&
        Object.prototype.hasOwnProperty.call(elem.attribs, name) &&
        elem.attribs[name] != null);
}
traversal.hasAttrib = hasAttrib;
/**
 * Get the tag name of an element.
 *
 * @param elem The element to get the name for.
 * @returns The tag name of `elem`.
 */
function getName(elem) {
    return elem.name;
}
traversal.getName = getName;
/**
 * Returns the next element sibling of a node.
 *
 * @param elem The element to get the next sibling of.
 * @returns `elem`'s next sibling that is a tag.
 */
function nextElementSibling(elem) {
    var _a;
    var next = elem.next;
    while (next !== null && !(0, domhandler_1$4.isTag)(next))
        (_a = next, next = _a.next);
    return next;
}
traversal.nextElementSibling = nextElementSibling;
/**
 * Returns the previous element sibling of a node.
 *
 * @param elem The element to get the previous sibling of.
 * @returns `elem`'s previous sibling that is a tag.
 */
function prevElementSibling(elem) {
    var _a;
    var prev = elem.prev;
    while (prev !== null && !(0, domhandler_1$4.isTag)(prev))
        (_a = prev, prev = _a.prev);
    return prev;
}
traversal.prevElementSibling = prevElementSibling;

var manipulation = {};

Object.defineProperty(manipulation, "__esModule", { value: true });
manipulation.prepend = manipulation.prependChild = manipulation.append = manipulation.appendChild = manipulation.replaceElement = manipulation.removeElement = void 0;
/**
 * Remove an element from the dom
 *
 * @param elem The element to be removed
 */
function removeElement(elem) {
    if (elem.prev)
        elem.prev.next = elem.next;
    if (elem.next)
        elem.next.prev = elem.prev;
    if (elem.parent) {
        var childs = elem.parent.children;
        childs.splice(childs.lastIndexOf(elem), 1);
    }
}
manipulation.removeElement = removeElement;
/**
 * Replace an element in the dom
 *
 * @param elem The element to be replaced
 * @param replacement The element to be added
 */
function replaceElement(elem, replacement) {
    var prev = (replacement.prev = elem.prev);
    if (prev) {
        prev.next = replacement;
    }
    var next = (replacement.next = elem.next);
    if (next) {
        next.prev = replacement;
    }
    var parent = (replacement.parent = elem.parent);
    if (parent) {
        var childs = parent.children;
        childs[childs.lastIndexOf(elem)] = replacement;
    }
}
manipulation.replaceElement = replaceElement;
/**
 * Append a child to an element.
 *
 * @param elem The element to append to.
 * @param child The element to be added as a child.
 */
function appendChild(elem, child) {
    removeElement(child);
    child.next = null;
    child.parent = elem;
    if (elem.children.push(child) > 1) {
        var sibling = elem.children[elem.children.length - 2];
        sibling.next = child;
        child.prev = sibling;
    }
    else {
        child.prev = null;
    }
}
manipulation.appendChild = appendChild;
/**
 * Append an element after another.
 *
 * @param elem The element to append after.
 * @param next The element be added.
 */
function append(elem, next) {
    removeElement(next);
    var parent = elem.parent;
    var currNext = elem.next;
    next.next = currNext;
    next.prev = elem;
    elem.next = next;
    next.parent = parent;
    if (currNext) {
        currNext.prev = next;
        if (parent) {
            var childs = parent.children;
            childs.splice(childs.lastIndexOf(currNext), 0, next);
        }
    }
    else if (parent) {
        parent.children.push(next);
    }
}
manipulation.append = append;
/**
 * Prepend a child to an element.
 *
 * @param elem The element to prepend before.
 * @param child The element to be added as a child.
 */
function prependChild(elem, child) {
    removeElement(child);
    child.parent = elem;
    child.prev = null;
    if (elem.children.unshift(child) !== 1) {
        var sibling = elem.children[1];
        sibling.prev = child;
        child.next = sibling;
    }
    else {
        child.next = null;
    }
}
manipulation.prependChild = prependChild;
/**
 * Prepend an element before another.
 *
 * @param elem The element to prepend before.
 * @param prev The element be added.
 */
function prepend(elem, prev) {
    removeElement(prev);
    var parent = elem.parent;
    if (parent) {
        var childs = parent.children;
        childs.splice(childs.indexOf(elem), 0, prev);
    }
    if (elem.prev) {
        elem.prev.next = prev;
    }
    prev.parent = parent;
    prev.prev = elem.prev;
    prev.next = elem;
    elem.prev = prev;
}
manipulation.prepend = prepend;

var querying = {};

Object.defineProperty(querying, "__esModule", { value: true });
querying.findAll = querying.existsOne = querying.findOne = querying.findOneChild = querying.find = querying.filter = void 0;
var domhandler_1$3 = lib$5;
/**
 * Search a node and its children for nodes passing a test function.
 *
 * @param test Function to test nodes on.
 * @param node Node to search. Will be included in the result set if it matches.
 * @param recurse Also consider child nodes.
 * @param limit Maximum number of nodes to return.
 * @returns All nodes passing `test`.
 */
function filter(test, node, recurse, limit) {
    if (recurse === void 0) { recurse = true; }
    if (limit === void 0) { limit = Infinity; }
    if (!Array.isArray(node))
        node = [node];
    return find(test, node, recurse, limit);
}
querying.filter = filter;
/**
 * Search an array of node and its children for nodes passing a test function.
 *
 * @param test Function to test nodes on.
 * @param nodes Array of nodes to search.
 * @param recurse Also consider child nodes.
 * @param limit Maximum number of nodes to return.
 * @returns All nodes passing `test`.
 */
function find(test, nodes, recurse, limit) {
    var result = [];
    for (var _i = 0, nodes_1 = nodes; _i < nodes_1.length; _i++) {
        var elem = nodes_1[_i];
        if (test(elem)) {
            result.push(elem);
            if (--limit <= 0)
                break;
        }
        if (recurse && (0, domhandler_1$3.hasChildren)(elem) && elem.children.length > 0) {
            var children = find(test, elem.children, recurse, limit);
            result.push.apply(result, children);
            limit -= children.length;
            if (limit <= 0)
                break;
        }
    }
    return result;
}
querying.find = find;
/**
 * Finds the first element inside of an array that matches a test function.
 *
 * @param test Function to test nodes on.
 * @param nodes Array of nodes to search.
 * @returns The first node in the array that passes `test`.
 */
function findOneChild(test, nodes) {
    return nodes.find(test);
}
querying.findOneChild = findOneChild;
/**
 * Finds one element in a tree that passes a test.
 *
 * @param test Function to test nodes on.
 * @param nodes Array of nodes to search.
 * @param recurse Also consider child nodes.
 * @returns The first child node that passes `test`.
 */
function findOne(test, nodes, recurse) {
    if (recurse === void 0) { recurse = true; }
    var elem = null;
    for (var i = 0; i < nodes.length && !elem; i++) {
        var checked = nodes[i];
        if (!(0, domhandler_1$3.isTag)(checked)) {
            continue;
        }
        else if (test(checked)) {
            elem = checked;
        }
        else if (recurse && checked.children.length > 0) {
            elem = findOne(test, checked.children);
        }
    }
    return elem;
}
querying.findOne = findOne;
/**
 * @param test Function to test nodes on.
 * @param nodes Array of nodes to search.
 * @returns Whether a tree of nodes contains at least one node passing a test.
 */
function existsOne(test, nodes) {
    return nodes.some(function (checked) {
        return (0, domhandler_1$3.isTag)(checked) &&
            (test(checked) ||
                (checked.children.length > 0 &&
                    existsOne(test, checked.children)));
    });
}
querying.existsOne = existsOne;
/**
 * Search and array of nodes and its children for nodes passing a test function.
 *
 * Same as `find`, only with less options, leading to reduced complexity.
 *
 * @param test Function to test nodes on.
 * @param nodes Array of nodes to search.
 * @returns All nodes passing `test`.
 */
function findAll(test, nodes) {
    var _a;
    var result = [];
    var stack = nodes.filter(domhandler_1$3.isTag);
    var elem;
    while ((elem = stack.shift())) {
        var children = (_a = elem.children) === null || _a === void 0 ? void 0 : _a.filter(domhandler_1$3.isTag);
        if (children && children.length > 0) {
            stack.unshift.apply(stack, children);
        }
        if (test(elem))
            result.push(elem);
    }
    return result;
}
querying.findAll = findAll;

var legacy = {};

Object.defineProperty(legacy, "__esModule", { value: true });
legacy.getElementsByTagType = legacy.getElementsByTagName = legacy.getElementById = legacy.getElements = legacy.testElement = void 0;
var domhandler_1$2 = lib$5;
var querying_1 = querying;
var Checks = {
    tag_name: function (name) {
        if (typeof name === "function") {
            return function (elem) { return (0, domhandler_1$2.isTag)(elem) && name(elem.name); };
        }
        else if (name === "*") {
            return domhandler_1$2.isTag;
        }
        return function (elem) { return (0, domhandler_1$2.isTag)(elem) && elem.name === name; };
    },
    tag_type: function (type) {
        if (typeof type === "function") {
            return function (elem) { return type(elem.type); };
        }
        return function (elem) { return elem.type === type; };
    },
    tag_contains: function (data) {
        if (typeof data === "function") {
            return function (elem) { return (0, domhandler_1$2.isText)(elem) && data(elem.data); };
        }
        return function (elem) { return (0, domhandler_1$2.isText)(elem) && elem.data === data; };
    },
};
/**
 * @param attrib Attribute to check.
 * @param value Attribute value to look for.
 * @returns A function to check whether the a node has an attribute with a particular value.
 */
function getAttribCheck(attrib, value) {
    if (typeof value === "function") {
        return function (elem) { return (0, domhandler_1$2.isTag)(elem) && value(elem.attribs[attrib]); };
    }
    return function (elem) { return (0, domhandler_1$2.isTag)(elem) && elem.attribs[attrib] === value; };
}
/**
 * @param a First function to combine.
 * @param b Second function to combine.
 * @returns A function taking a node and returning `true` if either
 * of the input functions returns `true` for the node.
 */
function combineFuncs(a, b) {
    return function (elem) { return a(elem) || b(elem); };
}
/**
 * @param options An object describing nodes to look for.
 * @returns A function executing all checks in `options` and returning `true`
 * if any of them match a node.
 */
function compileTest(options) {
    var funcs = Object.keys(options).map(function (key) {
        var value = options[key];
        return Object.prototype.hasOwnProperty.call(Checks, key)
            ? Checks[key](value)
            : getAttribCheck(key, value);
    });
    return funcs.length === 0 ? null : funcs.reduce(combineFuncs);
}
/**
 * @param options An object describing nodes to look for.
 * @param node The element to test.
 * @returns Whether the element matches the description in `options`.
 */
function testElement(options, node) {
    var test = compileTest(options);
    return test ? test(node) : true;
}
legacy.testElement = testElement;
/**
 * @param options An object describing nodes to look for.
 * @param nodes Nodes to search through.
 * @param recurse Also consider child nodes.
 * @param limit Maximum number of nodes to return.
 * @returns All nodes that match `options`.
 */
function getElements$1(options, nodes, recurse, limit) {
    if (limit === void 0) { limit = Infinity; }
    var test = compileTest(options);
    return test ? (0, querying_1.filter)(test, nodes, recurse, limit) : [];
}
legacy.getElements = getElements$1;
/**
 * @param id The unique ID attribute value to look for.
 * @param nodes Nodes to search through.
 * @param recurse Also consider child nodes.
 * @returns The node with the supplied ID.
 */
function getElementById(id, nodes, recurse) {
    if (recurse === void 0) { recurse = true; }
    if (!Array.isArray(nodes))
        nodes = [nodes];
    return (0, querying_1.findOne)(getAttribCheck("id", id), nodes, recurse);
}
legacy.getElementById = getElementById;
/**
 * @param tagName Tag name to search for.
 * @param nodes Nodes to search through.
 * @param recurse Also consider child nodes.
 * @param limit Maximum number of nodes to return.
 * @returns All nodes with the supplied `tagName`.
 */
function getElementsByTagName(tagName, nodes, recurse, limit) {
    if (recurse === void 0) { recurse = true; }
    if (limit === void 0) { limit = Infinity; }
    return (0, querying_1.filter)(Checks.tag_name(tagName), nodes, recurse, limit);
}
legacy.getElementsByTagName = getElementsByTagName;
/**
 * @param type Element type to look for.
 * @param nodes Nodes to search through.
 * @param recurse Also consider child nodes.
 * @param limit Maximum number of nodes to return.
 * @returns All nodes with the supplied `type`.
 */
function getElementsByTagType(type, nodes, recurse, limit) {
    if (recurse === void 0) { recurse = true; }
    if (limit === void 0) { limit = Infinity; }
    return (0, querying_1.filter)(Checks.tag_type(type), nodes, recurse, limit);
}
legacy.getElementsByTagType = getElementsByTagType;

var helpers = {};

Object.defineProperty(helpers, "__esModule", { value: true });
helpers.uniqueSort = helpers.compareDocumentPosition = helpers.removeSubsets = void 0;
var domhandler_1$1 = lib$5;
/**
 * Given an array of nodes, remove any member that is contained by another.
 *
 * @param nodes Nodes to filter.
 * @returns Remaining nodes that aren't subtrees of each other.
 */
function removeSubsets(nodes) {
    var idx = nodes.length;
    /*
     * Check if each node (or one of its ancestors) is already contained in the
     * array.
     */
    while (--idx >= 0) {
        var node = nodes[idx];
        /*
         * Remove the node if it is not unique.
         * We are going through the array from the end, so we only
         * have to check nodes that preceed the node under consideration in the array.
         */
        if (idx > 0 && nodes.lastIndexOf(node, idx - 1) >= 0) {
            nodes.splice(idx, 1);
            continue;
        }
        for (var ancestor = node.parent; ancestor; ancestor = ancestor.parent) {
            if (nodes.includes(ancestor)) {
                nodes.splice(idx, 1);
                break;
            }
        }
    }
    return nodes;
}
helpers.removeSubsets = removeSubsets;
/**
 * Compare the position of one node against another node in any other document.
 * The return value is a bitmask with the following values:
 *
 * Document order:
 * > There is an ordering, document order, defined on all the nodes in the
 * > document corresponding to the order in which the first character of the
 * > XML representation of each node occurs in the XML representation of the
 * > document after expansion of general entities. Thus, the document element
 * > node will be the first node. Element nodes occur before their children.
 * > Thus, document order orders element nodes in order of the occurrence of
 * > their start-tag in the XML (after expansion of entities). The attribute
 * > nodes of an element occur after the element and before its children. The
 * > relative order of attribute nodes is implementation-dependent./
 *
 * Source:
 * http://www.w3.org/TR/DOM-Level-3-Core/glossary.html#dt-document-order
 *
 * @param nodeA The first node to use in the comparison
 * @param nodeB The second node to use in the comparison
 * @returns A bitmask describing the input nodes' relative position.
 *
 * See http://dom.spec.whatwg.org/#dom-node-comparedocumentposition for
 * a description of these values.
 */
function compareDocumentPosition(nodeA, nodeB) {
    var aParents = [];
    var bParents = [];
    if (nodeA === nodeB) {
        return 0;
    }
    var current = (0, domhandler_1$1.hasChildren)(nodeA) ? nodeA : nodeA.parent;
    while (current) {
        aParents.unshift(current);
        current = current.parent;
    }
    current = (0, domhandler_1$1.hasChildren)(nodeB) ? nodeB : nodeB.parent;
    while (current) {
        bParents.unshift(current);
        current = current.parent;
    }
    var maxIdx = Math.min(aParents.length, bParents.length);
    var idx = 0;
    while (idx < maxIdx && aParents[idx] === bParents[idx]) {
        idx++;
    }
    if (idx === 0) {
        return 1 /* DISCONNECTED */;
    }
    var sharedParent = aParents[idx - 1];
    var siblings = sharedParent.children;
    var aSibling = aParents[idx];
    var bSibling = bParents[idx];
    if (siblings.indexOf(aSibling) > siblings.indexOf(bSibling)) {
        if (sharedParent === nodeB) {
            return 4 /* FOLLOWING */ | 16 /* CONTAINED_BY */;
        }
        return 4 /* FOLLOWING */;
    }
    if (sharedParent === nodeA) {
        return 2 /* PRECEDING */ | 8 /* CONTAINS */;
    }
    return 2 /* PRECEDING */;
}
helpers.compareDocumentPosition = compareDocumentPosition;
/**
 * Sort an array of nodes based on their relative position in the document and
 * remove any duplicate nodes. If the array contains nodes that do not belong
 * to the same document, sort order is unspecified.
 *
 * @param nodes Array of DOM nodes.
 * @returns Collection of unique nodes, sorted in document order.
 */
function uniqueSort(nodes) {
    nodes = nodes.filter(function (node, i, arr) { return !arr.includes(node, i + 1); });
    nodes.sort(function (a, b) {
        var relative = compareDocumentPosition(a, b);
        if (relative & 2 /* PRECEDING */) {
            return -1;
        }
        else if (relative & 4 /* FOLLOWING */) {
            return 1;
        }
        return 0;
    });
    return nodes;
}
helpers.uniqueSort = uniqueSort;

var feeds = {};

Object.defineProperty(feeds, "__esModule", { value: true });
feeds.getFeed = void 0;
var stringify_1 = stringify;
var legacy_1 = legacy;
/**
 * Get the feed object from the root of a DOM tree.
 *
 * @param doc - The DOM to to extract the feed from.
 * @returns The feed.
 */
function getFeed(doc) {
    var feedRoot = getOneElement$1(isValidFeed$1, doc);
    return !feedRoot
        ? null
        : feedRoot.name === "feed"
            ? getAtomFeed(feedRoot)
            : getRssFeed(feedRoot);
}
feeds.getFeed = getFeed;
/**
 * Parse an Atom feed.
 *
 * @param feedRoot The root of the feed.
 * @returns The parsed feed.
 */
function getAtomFeed(feedRoot) {
    var _a;
    var childs = feedRoot.children;
    var feed = {
        type: "atom",
        items: (0, legacy_1.getElementsByTagName)("entry", childs).map(function (item) {
            var _a;
            var children = item.children;
            var entry = { media: getMediaElements$1(children) };
            addConditionally$1(entry, "id", "id", children);
            addConditionally$1(entry, "title", "title", children);
            var href = (_a = getOneElement$1("link", children)) === null || _a === void 0 ? void 0 : _a.attribs.href;
            if (href) {
                entry.link = href;
            }
            var description = fetch$1("summary", children) || fetch$1("content", children);
            if (description) {
                entry.description = description;
            }
            var pubDate = fetch$1("updated", children);
            if (pubDate) {
                entry.pubDate = new Date(pubDate);
            }
            return entry;
        }),
    };
    addConditionally$1(feed, "id", "id", childs);
    addConditionally$1(feed, "title", "title", childs);
    var href = (_a = getOneElement$1("link", childs)) === null || _a === void 0 ? void 0 : _a.attribs.href;
    if (href) {
        feed.link = href;
    }
    addConditionally$1(feed, "description", "subtitle", childs);
    var updated = fetch$1("updated", childs);
    if (updated) {
        feed.updated = new Date(updated);
    }
    addConditionally$1(feed, "author", "email", childs, true);
    return feed;
}
/**
 * Parse a RSS feed.
 *
 * @param feedRoot The root of the feed.
 * @returns The parsed feed.
 */
function getRssFeed(feedRoot) {
    var _a, _b;
    var childs = (_b = (_a = getOneElement$1("channel", feedRoot.children)) === null || _a === void 0 ? void 0 : _a.children) !== null && _b !== void 0 ? _b : [];
    var feed = {
        type: feedRoot.name.substr(0, 3),
        id: "",
        items: (0, legacy_1.getElementsByTagName)("item", feedRoot.children).map(function (item) {
            var children = item.children;
            var entry = { media: getMediaElements$1(children) };
            addConditionally$1(entry, "id", "guid", children);
            addConditionally$1(entry, "title", "title", children);
            addConditionally$1(entry, "link", "link", children);
            addConditionally$1(entry, "description", "description", children);
            var pubDate = fetch$1("pubDate", children);
            if (pubDate)
                entry.pubDate = new Date(pubDate);
            return entry;
        }),
    };
    addConditionally$1(feed, "title", "title", childs);
    addConditionally$1(feed, "link", "link", childs);
    addConditionally$1(feed, "description", "description", childs);
    var updated = fetch$1("lastBuildDate", childs);
    if (updated) {
        feed.updated = new Date(updated);
    }
    addConditionally$1(feed, "author", "managingEditor", childs, true);
    return feed;
}
var MEDIA_KEYS_STRING = ["url", "type", "lang"];
var MEDIA_KEYS_INT = [
    "fileSize",
    "bitrate",
    "framerate",
    "samplingrate",
    "channels",
    "duration",
    "height",
    "width",
];
/**
 * Get all media elements of a feed item.
 *
 * @param where Nodes to search in.
 * @returns Media elements.
 */
function getMediaElements$1(where) {
    return (0, legacy_1.getElementsByTagName)("media:content", where).map(function (elem) {
        var attribs = elem.attribs;
        var media = {
            medium: attribs.medium,
            isDefault: !!attribs.isDefault,
        };
        for (var _i = 0, MEDIA_KEYS_STRING_1 = MEDIA_KEYS_STRING; _i < MEDIA_KEYS_STRING_1.length; _i++) {
            var attrib = MEDIA_KEYS_STRING_1[_i];
            if (attribs[attrib]) {
                media[attrib] = attribs[attrib];
            }
        }
        for (var _a = 0, MEDIA_KEYS_INT_1 = MEDIA_KEYS_INT; _a < MEDIA_KEYS_INT_1.length; _a++) {
            var attrib = MEDIA_KEYS_INT_1[_a];
            if (attribs[attrib]) {
                media[attrib] = parseInt(attribs[attrib], 10);
            }
        }
        if (attribs.expression) {
            media.expression =
                attribs.expression;
        }
        return media;
    });
}
/**
 * Get one element by tag name.
 *
 * @param tagName Tag name to look for
 * @param node Node to search in
 * @returns The element or null
 */
function getOneElement$1(tagName, node) {
    return (0, legacy_1.getElementsByTagName)(tagName, node, true, 1)[0];
}
/**
 * Get the text content of an element with a certain tag name.
 *
 * @param tagName Tag name to look for.
 * @param where  Node to search in.
 * @param recurse Whether to recurse into child nodes.
 * @returns The text content of the element.
 */
function fetch$1(tagName, where, recurse) {
    if (recurse === void 0) { recurse = false; }
    return (0, stringify_1.textContent)((0, legacy_1.getElementsByTagName)(tagName, where, recurse, 1)).trim();
}
/**
 * Adds a property to an object if it has a value.
 *
 * @param obj Object to be extended
 * @param prop Property name
 * @param tagName Tag name that contains the conditionally added property
 * @param where Element to search for the property
 * @param recurse Whether to recurse into child nodes.
 */
function addConditionally$1(obj, prop, tagName, where, recurse) {
    if (recurse === void 0) { recurse = false; }
    var val = fetch$1(tagName, where, recurse);
    if (val)
        obj[prop] = val;
}
/**
 * Checks if an element is a feed root node.
 *
 * @param value The name of the element to check.
 * @returns Whether an element is a feed root node.
 */
function isValidFeed$1(value) {
    return value === "rss" || value === "feed" || value === "rdf:RDF";
}

(function (exports) {
	var __createBinding = (commonjsGlobal && commonjsGlobal.__createBinding) || (Object.create ? (function(o, m, k, k2) {
	    if (k2 === undefined) k2 = k;
	    Object.defineProperty(o, k2, { enumerable: true, get: function() { return m[k]; } });
	}) : (function(o, m, k, k2) {
	    if (k2 === undefined) k2 = k;
	    o[k2] = m[k];
	}));
	var __exportStar = (commonjsGlobal && commonjsGlobal.__exportStar) || function(m, exports) {
	    for (var p in m) if (p !== "default" && !Object.prototype.hasOwnProperty.call(exports, p)) __createBinding(exports, m, p);
	};
	Object.defineProperty(exports, "__esModule", { value: true });
	exports.hasChildren = exports.isDocument = exports.isComment = exports.isText = exports.isCDATA = exports.isTag = void 0;
	__exportStar(stringify, exports);
	__exportStar(traversal, exports);
	__exportStar(manipulation, exports);
	__exportStar(querying, exports);
	__exportStar(legacy, exports);
	__exportStar(helpers, exports);
	__exportStar(feeds, exports);
	/** @deprecated Use these methods from `domhandler` directly. */
	var domhandler_1 = lib$5;
	Object.defineProperty(exports, "isTag", { enumerable: true, get: function () { return domhandler_1.isTag; } });
	Object.defineProperty(exports, "isCDATA", { enumerable: true, get: function () { return domhandler_1.isCDATA; } });
	Object.defineProperty(exports, "isText", { enumerable: true, get: function () { return domhandler_1.isText; } });
	Object.defineProperty(exports, "isComment", { enumerable: true, get: function () { return domhandler_1.isComment; } });
	Object.defineProperty(exports, "isDocument", { enumerable: true, get: function () { return domhandler_1.isDocument; } });
	Object.defineProperty(exports, "hasChildren", { enumerable: true, get: function () { return domhandler_1.hasChildren; } }); 
} (lib$2));

var __extends = (commonjsGlobal && commonjsGlobal.__extends) || (function () {
    var extendStatics = function (d, b) {
        extendStatics = Object.setPrototypeOf ||
            ({ __proto__: [] } instanceof Array && function (d, b) { d.__proto__ = b; }) ||
            function (d, b) { for (var p in b) if (Object.prototype.hasOwnProperty.call(b, p)) d[p] = b[p]; };
        return extendStatics(d, b);
    };
    return function (d, b) {
        if (typeof b !== "function" && b !== null)
            throw new TypeError("Class extends value " + String(b) + " is not a constructor or null");
        extendStatics(d, b);
        function __() { this.constructor = d; }
        d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
    };
})();
var __createBinding = (commonjsGlobal && commonjsGlobal.__createBinding) || (Object.create ? (function(o, m, k, k2) {
    if (k2 === undefined) k2 = k;
    Object.defineProperty(o, k2, { enumerable: true, get: function() { return m[k]; } });
}) : (function(o, m, k, k2) {
    if (k2 === undefined) k2 = k;
    o[k2] = m[k];
}));
var __setModuleDefault = (commonjsGlobal && commonjsGlobal.__setModuleDefault) || (Object.create ? (function(o, v) {
    Object.defineProperty(o, "default", { enumerable: true, value: v });
}) : function(o, v) {
    o["default"] = v;
});
var __importStar = (commonjsGlobal && commonjsGlobal.__importStar) || function (mod) {
    if (mod && mod.__esModule) return mod;
    var result = {};
    if (mod != null) for (var k in mod) if (k !== "default" && Object.prototype.hasOwnProperty.call(mod, k)) __createBinding(result, mod, k);
    __setModuleDefault(result, mod);
    return result;
};
var __importDefault = (commonjsGlobal && commonjsGlobal.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(FeedHandler$1, "__esModule", { value: true });
FeedHandler$1.parseFeed = FeedHandler$1.FeedHandler = void 0;
var domhandler_1 = __importDefault(lib$5);
var DomUtils = __importStar(lib$2);
var Parser_1 = Parser$1;
var FeedItemMediaMedium;
(function (FeedItemMediaMedium) {
    FeedItemMediaMedium[FeedItemMediaMedium["image"] = 0] = "image";
    FeedItemMediaMedium[FeedItemMediaMedium["audio"] = 1] = "audio";
    FeedItemMediaMedium[FeedItemMediaMedium["video"] = 2] = "video";
    FeedItemMediaMedium[FeedItemMediaMedium["document"] = 3] = "document";
    FeedItemMediaMedium[FeedItemMediaMedium["executable"] = 4] = "executable";
})(FeedItemMediaMedium || (FeedItemMediaMedium = {}));
var FeedItemMediaExpression;
(function (FeedItemMediaExpression) {
    FeedItemMediaExpression[FeedItemMediaExpression["sample"] = 0] = "sample";
    FeedItemMediaExpression[FeedItemMediaExpression["full"] = 1] = "full";
    FeedItemMediaExpression[FeedItemMediaExpression["nonstop"] = 2] = "nonstop";
})(FeedItemMediaExpression || (FeedItemMediaExpression = {}));
// TODO: Consume data as it is coming in
var FeedHandler = /** @class */ (function (_super) {
    __extends(FeedHandler, _super);
    /**
     *
     * @param callback
     * @param options
     */
    function FeedHandler(callback, options) {
        var _this = this;
        if (typeof callback === "object") {
            callback = undefined;
            options = callback;
        }
        _this = _super.call(this, callback, options) || this;
        return _this;
    }
    FeedHandler.prototype.onend = function () {
        var _a, _b;
        var feedRoot = getOneElement(isValidFeed, this.dom);
        if (!feedRoot) {
            this.handleCallback(new Error("couldn't find root of feed"));
            return;
        }
        var feed = {};
        if (feedRoot.name === "feed") {
            var childs = feedRoot.children;
            feed.type = "atom";
            addConditionally(feed, "id", "id", childs);
            addConditionally(feed, "title", "title", childs);
            var href = getAttribute("href", getOneElement("link", childs));
            if (href) {
                feed.link = href;
            }
            addConditionally(feed, "description", "subtitle", childs);
            var updated = fetch("updated", childs);
            if (updated) {
                feed.updated = new Date(updated);
            }
            addConditionally(feed, "author", "email", childs, true);
            feed.items = getElements("entry", childs).map(function (item) {
                var entry = {};
                var children = item.children;
                addConditionally(entry, "id", "id", children);
                addConditionally(entry, "title", "title", children);
                var href = getAttribute("href", getOneElement("link", children));
                if (href) {
                    entry.link = href;
                }
                var description = fetch("summary", children) || fetch("content", children);
                if (description) {
                    entry.description = description;
                }
                var pubDate = fetch("updated", children);
                if (pubDate) {
                    entry.pubDate = new Date(pubDate);
                }
                entry.media = getMediaElements(children);
                return entry;
            });
        }
        else {
            var childs = (_b = (_a = getOneElement("channel", feedRoot.children)) === null || _a === void 0 ? void 0 : _a.children) !== null && _b !== void 0 ? _b : [];
            feed.type = feedRoot.name.substr(0, 3);
            feed.id = "";
            addConditionally(feed, "title", "title", childs);
            addConditionally(feed, "link", "link", childs);
            addConditionally(feed, "description", "description", childs);
            var updated = fetch("lastBuildDate", childs);
            if (updated) {
                feed.updated = new Date(updated);
            }
            addConditionally(feed, "author", "managingEditor", childs, true);
            feed.items = getElements("item", feedRoot.children).map(function (item) {
                var entry = {};
                var children = item.children;
                addConditionally(entry, "id", "guid", children);
                addConditionally(entry, "title", "title", children);
                addConditionally(entry, "link", "link", children);
                addConditionally(entry, "description", "description", children);
                var pubDate = fetch("pubDate", children);
                if (pubDate)
                    entry.pubDate = new Date(pubDate);
                entry.media = getMediaElements(children);
                return entry;
            });
        }
        this.feed = feed;
        this.handleCallback(null);
    };
    return FeedHandler;
}(domhandler_1.default));
FeedHandler$1.FeedHandler = FeedHandler;
function getMediaElements(where) {
    return getElements("media:content", where).map(function (elem) {
        var media = {
            medium: elem.attribs.medium,
            isDefault: !!elem.attribs.isDefault,
        };
        if (elem.attribs.url) {
            media.url = elem.attribs.url;
        }
        if (elem.attribs.fileSize) {
            media.fileSize = parseInt(elem.attribs.fileSize, 10);
        }
        if (elem.attribs.type) {
            media.type = elem.attribs.type;
        }
        if (elem.attribs.expression) {
            media.expression = elem.attribs
                .expression;
        }
        if (elem.attribs.bitrate) {
            media.bitrate = parseInt(elem.attribs.bitrate, 10);
        }
        if (elem.attribs.framerate) {
            media.framerate = parseInt(elem.attribs.framerate, 10);
        }
        if (elem.attribs.samplingrate) {
            media.samplingrate = parseInt(elem.attribs.samplingrate, 10);
        }
        if (elem.attribs.channels) {
            media.channels = parseInt(elem.attribs.channels, 10);
        }
        if (elem.attribs.duration) {
            media.duration = parseInt(elem.attribs.duration, 10);
        }
        if (elem.attribs.height) {
            media.height = parseInt(elem.attribs.height, 10);
        }
        if (elem.attribs.width) {
            media.width = parseInt(elem.attribs.width, 10);
        }
        if (elem.attribs.lang) {
            media.lang = elem.attribs.lang;
        }
        return media;
    });
}
function getElements(tagName, where) {
    return DomUtils.getElementsByTagName(tagName, where, true);
}
function getOneElement(tagName, node) {
    return DomUtils.getElementsByTagName(tagName, node, true, 1)[0];
}
function fetch(tagName, where, recurse) {
    if (recurse === void 0) { recurse = false; }
    return DomUtils.getText(DomUtils.getElementsByTagName(tagName, where, recurse, 1)).trim();
}
function getAttribute(name, elem) {
    if (!elem) {
        return null;
    }
    var attribs = elem.attribs;
    return attribs[name];
}
function addConditionally(obj, prop, what, where, recurse) {
    if (recurse === void 0) { recurse = false; }
    var tmp = fetch(what, where, recurse);
    if (tmp)
        obj[prop] = tmp;
}
function isValidFeed(value) {
    return value === "rss" || value === "feed" || value === "rdf:RDF";
}
/**
 * Parse a feed.
 *
 * @param feed The feed that should be parsed, as a string.
 * @param options Optionally, options for parsing. When using this option, you should set `xmlMode` to `true`.
 */
function parseFeed(feed, options) {
    if (options === void 0) { options = { xmlMode: true }; }
    var handler = new FeedHandler(options);
    new Parser_1.Parser(handler, options).end(feed);
    return handler.feed;
}
FeedHandler$1.parseFeed = parseFeed;

(function (exports) {
	var __createBinding = (commonjsGlobal && commonjsGlobal.__createBinding) || (Object.create ? (function(o, m, k, k2) {
	    if (k2 === undefined) k2 = k;
	    Object.defineProperty(o, k2, { enumerable: true, get: function() { return m[k]; } });
	}) : (function(o, m, k, k2) {
	    if (k2 === undefined) k2 = k;
	    o[k2] = m[k];
	}));
	var __setModuleDefault = (commonjsGlobal && commonjsGlobal.__setModuleDefault) || (Object.create ? (function(o, v) {
	    Object.defineProperty(o, "default", { enumerable: true, value: v });
	}) : function(o, v) {
	    o["default"] = v;
	});
	var __importStar = (commonjsGlobal && commonjsGlobal.__importStar) || function (mod) {
	    if (mod && mod.__esModule) return mod;
	    var result = {};
	    if (mod != null) for (var k in mod) if (k !== "default" && Object.prototype.hasOwnProperty.call(mod, k)) __createBinding(result, mod, k);
	    __setModuleDefault(result, mod);
	    return result;
	};
	var __exportStar = (commonjsGlobal && commonjsGlobal.__exportStar) || function(m, exports) {
	    for (var p in m) if (p !== "default" && !Object.prototype.hasOwnProperty.call(exports, p)) __createBinding(exports, m, p);
	};
	var __importDefault = (commonjsGlobal && commonjsGlobal.__importDefault) || function (mod) {
	    return (mod && mod.__esModule) ? mod : { "default": mod };
	};
	Object.defineProperty(exports, "__esModule", { value: true });
	exports.RssHandler = exports.DefaultHandler = exports.DomUtils = exports.ElementType = exports.Tokenizer = exports.createDomStream = exports.parseDOM = exports.parseDocument = exports.DomHandler = exports.Parser = void 0;
	var Parser_1 = Parser$1;
	Object.defineProperty(exports, "Parser", { enumerable: true, get: function () { return Parser_1.Parser; } });
	var domhandler_1 = lib$5;
	Object.defineProperty(exports, "DomHandler", { enumerable: true, get: function () { return domhandler_1.DomHandler; } });
	Object.defineProperty(exports, "DefaultHandler", { enumerable: true, get: function () { return domhandler_1.DomHandler; } });
	// Helper methods
	/**
	 * Parses the data, returns the resulting document.
	 *
	 * @param data The data that should be parsed.
	 * @param options Optional options for the parser and DOM builder.
	 */
	function parseDocument(data, options) {
	    var handler = new domhandler_1.DomHandler(undefined, options);
	    new Parser_1.Parser(handler, options).end(data);
	    return handler.root;
	}
	exports.parseDocument = parseDocument;
	/**
	 * Parses data, returns an array of the root nodes.
	 *
	 * Note that the root nodes still have a `Document` node as their parent.
	 * Use `parseDocument` to get the `Document` node instead.
	 *
	 * @param data The data that should be parsed.
	 * @param options Optional options for the parser and DOM builder.
	 * @deprecated Use `parseDocument` instead.
	 */
	function parseDOM(data, options) {
	    return parseDocument(data, options).children;
	}
	exports.parseDOM = parseDOM;
	/**
	 * Creates a parser instance, with an attached DOM handler.
	 *
	 * @param cb A callback that will be called once parsing has been completed.
	 * @param options Optional options for the parser and DOM builder.
	 * @param elementCb An optional callback that will be called every time a tag has been completed inside of the DOM.
	 */
	function createDomStream(cb, options, elementCb) {
	    var handler = new domhandler_1.DomHandler(cb, options, elementCb);
	    return new Parser_1.Parser(handler, options);
	}
	exports.createDomStream = createDomStream;
	var Tokenizer_1 = Tokenizer$1;
	Object.defineProperty(exports, "Tokenizer", { enumerable: true, get: function () { return __importDefault(Tokenizer_1).default; } });
	var ElementType = __importStar(lib$4);
	exports.ElementType = ElementType;
	/*
	 * All of the following exports exist for backwards-compatibility.
	 * They should probably be removed eventually.
	 */
	__exportStar(FeedHandler$1, exports);
	exports.DomUtils = __importStar(lib$2);
	var FeedHandler_1 = FeedHandler$1;
	Object.defineProperty(exports, "RssHandler", { enumerable: true, get: function () { return FeedHandler_1.FeedHandler; } }); 
} (lib$3));

const merge$1 = cjs;

/**
 * Given a list of class and ID selectors (prefixed with '.' and '#'),
 * return them as separate lists of names without prefixes.
 *
 * @param { string[] } selectors Class and ID selectors (`[".class", "#id"]` etc).
 * @returns { { classes: string[], ids: string[] } }
 */
function splitClassesAndIds$1 (selectors) {
  const classes = [];
  const ids = [];
  for (const selector of selectors) {
    if (selector.startsWith('.')) {
      classes.push(selector.substring(1));
    } else if (selector.startsWith('#')) {
      ids.push(selector.substring(1));
    }
  }
  return { classes: classes, ids: ids };
}

/**
 * Make a recursive function that will only run to a given depth
 * and switches to an alternative function at that depth. \
 * No limitation if `n` is `undefined` (Just wraps `f` in that case).
 *
 * @param   { number | undefined } n   Allowed depth of recursion. `undefined` for no limitation.
 * @param   { Function }           f   Function that accepts recursive callback as the first argument.
 * @param   { Function }           [g] Function to run instead, when maximum depth was reached. Do nothing by default.
 * @returns { Function }
 */
function limitedDepthRecursive$1 (n, f, g = () => undefined) {
  if (n === undefined) {
    const f1 = function (...args) { return f(f1, ...args); };
    return f1;
  }
  if (n >= 0) {
    return function (...args) { return f(limitedDepthRecursive$1(n - 1, f, g), ...args); };
  }
  return g;
}

/**
 * Convert a number into alphabetic sequence representation (Sequence without zeroes).
 *
 * For example: `a, ..., z, aa, ..., zz, aaa, ...`.
 *
 * @param   { number } num              Number to convert. Must be >= 1.
 * @param   { string } [baseChar = 'a'] Character for 1 in the sequence.
 * @param   { number } [base = 26]      Number of characters in the sequence.
 * @returns { string }
 */
function numberToLetterSequence$1 (num, baseChar = 'a', base = 26) {
  const digits = [];
  do {
    num -= 1;
    digits.push(num % base);
    num = (num / base) >> 0; // quick `floor`
  } while (num > 0);
  const baseCode = baseChar.charCodeAt(0);
  return digits
    .reverse()
    .map(n => String.fromCharCode(baseCode + n))
    .join('');
}

const I = ['I', 'X', 'C', 'M'];
const V = ['V', 'L', 'D'];

/**
 * Convert a number to it's Roman representation. No large numbers extension.
 *
 * @param   { number } num Number to convert. `0 < num <= 3999`.
 * @returns { string }
 */
function numberToRoman$1 (num) {
  return [...(num) + '']
    .map(n => +n)
    .reverse()
    .map((v, i) => ((v % 5 < 4)
      ? (v < 5 ? '' : V[i]) + I[i].repeat(v % 5)
      : I[i] + (v < 5 ? V[i] : I[i + 1])))
    .reverse()
    .join('');
}

/**
 * Return the same string or a substring with the given character occurences removed from each end if any.
 *
 * @param   { string } str  A string to trim.
 * @param   { string } char A character to be trimmed.
 * @returns { string }
 */
function trimCharacter$2 (str, char) {
  let start = 0;
  let end = str.length;
  while (start < end && str[start] === char) { ++start; }
  while (end > start && str[end - 1] === char) { --end; }
  return (start > 0 || end < str.length)
    ? str.substring(start, end)
    : str;
}

/**
 * Get a nested property from an object.
 *
 * @param   { object }   obj  The object to query for the value.
 * @param   { string[] } path The path to the property.
 * @returns { any }
 */
function get$2 (obj, path) {
  for (const key of path) {
    if (!obj) { return undefined; }
    obj = obj[key];
  }
  return obj;
}

/**
 * Deduplicate an array by a given key callback.
 * Item properties are merged recursively and with the preference for last defined values.
 * Of items with the same key, merged item takes the place of the last item,
 * others are omitted.
 *
 * @param { any[] } items An array to deduplicate.
 * @param { (x: any) => string } getKey Callback to get a value that distinguishes unique items.
 * @returns { any[] }
 */
function mergeDuplicatesPreferLast$1 (items, getKey) {
  const map = new Map();
  for (let i = items.length; i-- > 0;) {
    const item = items[i];
    const key = getKey(item);
    map.set(
      key,
      (map.has(key))
        ? merge$1(item, map.get(key), { arrayMerge: overwriteMerge$1 })
        : item
    );
  }
  return [...map.values()].reverse();
}

const overwriteMerge$1 = (acc, src, options) => [...src];

var helper = {
  get: get$2,
  limitedDepthRecursive: limitedDepthRecursive$1,
  mergeDuplicatesPreferLast: mergeDuplicatesPreferLast$1,
  numberToLetterSequence: numberToLetterSequence$1,
  numberToRoman: numberToRoman$1,
  splitClassesAndIds: splitClassesAndIds$1,
  trimCharacter: trimCharacter$2
};

// eslint-disable-next-line import/no-unassigned-import


/**
 * Helps to build text from words.
 */
let InlineTextBuilder$1 = class InlineTextBuilder {
  /**
   * Creates an instance of InlineTextBuilder.
   *
   * If `maxLineLength` is not provided then it is either `options.wordwrap` or unlimited.
   *
   * @param { Options } options           HtmlToText options.
   * @param { number }  [ maxLineLength ] This builder will try to wrap text to fit this line length.
   */
  constructor (options, maxLineLength = undefined) {
    /** @type { string[][] } */
    this.lines = [];
    /** @type { string[] }   */
    this.nextLineWords = [];
    this.maxLineLength = maxLineLength || options.wordwrap || Number.MAX_VALUE;
    this.nextLineAvailableChars = this.maxLineLength;
    this.wrapCharacters = options.longWordSplit.wrapCharacters || [];
    this.forceWrapOnLimit = options.longWordSplit.forceWrapOnLimit || false;

    this.stashedSpace = false;
    this.wordBreakOpportunity = false;
  }

  /**
   * Add a new word.
   *
   * @param { string } word A word to add.
   */
  pushWord (word) {
    if (this.nextLineAvailableChars <= 0) {
      this.startNewLine();
    }
    const isLineStart = this.nextLineWords.length === 0;
    const cost = word.length + (isLineStart ? 0 : 1);
    if (cost <= this.nextLineAvailableChars) { // Fits into available budget

      this.nextLineWords.push(word);
      this.nextLineAvailableChars -= cost;

    } else { // Does not fit - try to split the word

      // The word is moved to a new line - prefer to wrap between words.
      const [first, ...rest] = this.splitLongWord(word);
      if (!isLineStart) { this.startNewLine(); }
      this.nextLineWords.push(first);
      this.nextLineAvailableChars -= first.length;
      for (const part of rest) {
        this.startNewLine();
        this.nextLineWords.push(part);
        this.nextLineAvailableChars -= part.length;
      }

    }
  }

  /**
   * Pop a word from the currently built line.
   * This doesn't affect completed lines.
   *
   * @returns { string }
   */
  popWord () {
    const lastWord = this.nextLineWords.pop();
    if (lastWord !== undefined) {
      const isLineStart = this.nextLineWords.length === 0;
      const cost = lastWord.length + (isLineStart ? 0 : 1);
      this.nextLineAvailableChars += cost;
    }
    return lastWord;
  }

  /**
   * Concat a word to the last word already in the builder.
   * Adds a new word in case there are no words yet in the last line.
   *
   * @param { string } word A word to be concatenated.
   */
  concatWord (word) {
    if (this.wordBreakOpportunity && word.length > this.nextLineAvailableChars) {
      this.pushWord(word);
      this.wordBreakOpportunity = false;
    } else {
      const lastWord = this.popWord();
      this.pushWord((lastWord) ? lastWord.concat(word) : word);
    }
  }

  /**
   * Add current line (and more empty lines if provided argument > 1) to the list of complete lines and start a new one.
   *
   * @param { number } n Number of line breaks that will be added to the resulting string.
   */
  startNewLine (n = 1) {
    this.lines.push(this.nextLineWords);
    if (n > 1) {
      this.lines.push(...Array.from({ length: n - 1 }, () => []));
    }
    this.nextLineWords = [];
    this.nextLineAvailableChars = this.maxLineLength;
  }

  /**
   * No words in this builder.
   *
   * @returns { boolean }
   */
  isEmpty () {
    return this.lines.length === 0
        && this.nextLineWords.length === 0;
  }

  clear () {
    this.lines.length = 0;
    this.nextLineWords.length = 0;
    this.nextLineAvailableChars = this.maxLineLength;
  }

  /**
   * Join all lines of words inside the InlineTextBuilder into a complete string.
   *
   * @returns { string }
   */
  toString () {
    return [...this.lines, this.nextLineWords]
      .map(words => words.join(' '))
      .join('\n');
  }

  /**
   * Split a long word up to fit within the word wrap limit.
   * Use either a character to split looking back from the word wrap limit,
   * or truncate to the word wrap limit.
   *
   * @param   { string }   word Input word.
   * @returns { string[] }      Parts of the word.
   */
  splitLongWord (word) {
    const parts = [];
    let idx = 0;
    while (word.length > this.maxLineLength) {

      const firstLine = word.substring(0, this.maxLineLength);
      const remainingChars = word.substring(this.maxLineLength);

      const splitIndex = firstLine.lastIndexOf(this.wrapCharacters[idx]);

      if (splitIndex > -1) { // Found a character to split on

        word = firstLine.substring(splitIndex + 1) + remainingChars;
        parts.push(firstLine.substring(0, splitIndex + 1));

      } else { // Not found a character to split on

        idx++;
        if (idx < this.wrapCharacters.length) { // There is next character to try

          word = firstLine + remainingChars;

        } else { // No more characters to try

          if (this.forceWrapOnLimit) {
            parts.push(firstLine);
            word = remainingChars;
            if (word.length > this.maxLineLength) {
              continue;
            }
          } else {
            word = firstLine + remainingChars;
          }
          break;

        }

      }

    }
    parts.push(word); // Add remaining part to array
    return parts;
  }
};

var inlineTextBuilder = { InlineTextBuilder: InlineTextBuilder$1 };

/* eslint-disable max-classes-per-file */

const { InlineTextBuilder } = inlineTextBuilder;


let StackItem$1 = class StackItem {
  constructor (next = null) { this.next = next; }

  getRoot () { return (this.next) ? this.next : this; }
};

let BlockStackItem$1 = class BlockStackItem extends StackItem$1 {
  constructor (options, next = null, leadingLineBreaks = 1, maxLineLength = undefined) {
    super(next);
    this.leadingLineBreaks = leadingLineBreaks;
    this.inlineTextBuilder = new InlineTextBuilder(options, maxLineLength);
    this.rawText = '';
    this.stashedLineBreaks = 0;
    this.isPre = next && next.isPre;
  }
};

let TableStackItem$1 = class TableStackItem extends StackItem$1 {
  constructor (next = null) {
    super(next);
    this.rows = [];
    this.isPre = next && next.isPre;
  }
};

let TableRowStackItem$1 = class TableRowStackItem extends StackItem$1 {
  constructor (next = null) {
    super(next);
    this.cells = [];
    this.isPre = next && next.isPre;
  }
};

let TableCellStackItem$1 = class TableCellStackItem extends StackItem$1 {
  constructor (options, next = null, maxColumnWidth = undefined) {
    super(next);
    this.inlineTextBuilder = new InlineTextBuilder(options, maxColumnWidth);
    this.rawText = '';
    this.stashedLineBreaks = 0;
    this.isPre = next && next.isPre;
  }
};

let TransformerStackItem$1 = class TransformerStackItem extends StackItem$1 {
  constructor (next = null, transform) {
    super(next);
    this.transform = transform;
  }
};

var stackItem = {
  BlockStackItem: BlockStackItem$1,
  StackItem: StackItem$1,
  TableCellStackItem: TableCellStackItem$1,
  TableRowStackItem: TableRowStackItem$1,
  TableStackItem: TableStackItem$1,
  TransformerStackItem: TransformerStackItem$1,
};

function getRow (matrix, j) {
  if (!matrix[j]) { matrix[j] = []; }
  return matrix[j];
}

function findFirstVacantIndex (row, x = 0) {
  while (row[x]) { x++; }
  return x;
}

function transposeInPlace (matrix, maxSize) {
  for (let i = 0; i < maxSize; i++) {
    const rowI = getRow(matrix, i);
    for (let j = 0; j < i; j++) {
      const rowJ = getRow(matrix, j);
      const temp = rowI[j];
      rowI[j] = rowJ[i];
      rowJ[i] = temp;
    }
  }
}

function putCellIntoLayout (cell, layout, baseRow, baseCol) {
  for (let r = 0; r < cell.rowspan; r++) {
    const layoutRow = getRow(layout, baseRow + r);
    for (let c = 0; c < cell.colspan; c++) {
      layoutRow[baseCol + c] = cell;
    }
  }
}

function updateOffset (offsets, base, span, value) {
  offsets[base + span] = Math.max(
    offsets[base + span] || 0,
    offsets[base] + value
  );
}

/**
 * @typedef { object } TablePrinterCell
 * Cell definition for the table printer.
 *
 * @property { number } colspan Number of columns this cell occupies.
 * @property { number } rowspan Number of rows this cell occupies.
 * @property { string } text Cell contents (pre-wrapped).
 */

/**
 * Render a table into string.
 * Cells can contain multiline text and span across multiple rows and columns.
 *
 * Modifies cells to add lines array.
 *
 * @param { TablePrinterCell[][] } tableRows Table to render.
 * @param { number } rowSpacing Number of spaces between columns.
 * @param { number } colSpacing Number of empty lines between rows.
 * @returns { string }
 */
function tableToString$1 (tableRows, rowSpacing, colSpacing) {
  const layout = [];
  let colNumber = 0;
  const rowNumber = tableRows.length;
  const rowOffsets = [0];
  // Fill the layout table and row offsets row-by-row.
  for (let j = 0; j < rowNumber; j++) {
    const layoutRow = getRow(layout, j);
    const cells = tableRows[j];
    let x = 0;
    for (let i = 0; i < cells.length; i++) {
      const cell = cells[i];
      x = findFirstVacantIndex(layoutRow, x);
      putCellIntoLayout(cell, layout, j, x);
      x += cell.colspan;
      cell.lines = cell.text.split('\n');
      const cellHeight = cell.lines.length;
      updateOffset(rowOffsets, j, cell.rowspan, cellHeight + rowSpacing);
    }
    colNumber = (layoutRow.length > colNumber) ? layoutRow.length : colNumber;
  }

  transposeInPlace(layout, (rowNumber > colNumber) ? rowNumber : colNumber);

  const outputLines = [];
  const colOffsets = [0];
  // Fill column offsets and output lines column-by-column.
  for (let x = 0; x < colNumber; x++) {
    let y = 0;
    let cell;
    while (y < rowNumber && (cell = layout[x][y])) {
      if (!cell.rendered) {
        let cellWidth = 0;
        for (let j = 0; j < cell.lines.length; j++) {
          const line = cell.lines[j];
          const lineOffset = rowOffsets[y] + j;
          outputLines[lineOffset] = (outputLines[lineOffset] || '').padEnd(colOffsets[x]) + line;
          cellWidth = (line.length > cellWidth) ? line.length : cellWidth;
        }
        updateOffset(colOffsets, x, cell.colspan, cellWidth + colSpacing);
        cell.rendered = true;
      }
      y += cell.rowspan;
    }
  }

  return outputLines.join('\n');
}

var tablePrinter = { tableToString: tableToString$1 };

// eslint-disable-next-line import/no-unassigned-import



function charactersToCodes (str) {
  return [...str]
    .map(c => '\\u' + c.charCodeAt(0).toString(16).padStart(4, '0'))
    .join('');
}

/**
 * Helps to handle HTML whitespaces.
 *
 * @class WhitespaceProcessor
 */
let WhitespaceProcessor$1 = class WhitespaceProcessor {

  /**
   * Creates an instance of WhitespaceProcessor.
   *
   * @param { Options } options    HtmlToText options.
   * @memberof WhitespaceProcessor
   */
  constructor (options) {
    this.whitespaceChars = (options.preserveNewlines)
      ? options.whitespaceCharacters.replace(/\n/g, '')
      : options.whitespaceCharacters;
    const whitespaceCodes = charactersToCodes(this.whitespaceChars);
    this.leadingWhitespaceRe = new RegExp(`^[${whitespaceCodes}]`);
    this.trailingWhitespaceRe = new RegExp(`[${whitespaceCodes}]$`);
    this.allWhitespaceOrEmptyRe = new RegExp(`^[${whitespaceCodes}]*$`);
    this.newlineOrNonWhitespaceRe = new RegExp(`(\\n|[^\\n${whitespaceCodes}])`, 'g');

    if (options.preserveNewlines) {

      const wordOrNewlineRe = new RegExp(`\\n|[^\\n${whitespaceCodes}]+`, 'gm');

      /**
       * Shrink whitespaces and wrap text, add to the builder.
       *
       * @param { string }                  text              Input text.
       * @param { InlineTextBuilder }       inlineTextBuilder A builder to receive processed text.
       * @param { (str: string) => string } [ transform ]     A transform to be applied to words.
       */
      this.shrinkWrapAdd = function (text, inlineTextBuilder, transform = (str => str)) {
        if (!text) { return; }
        const previouslyStashedSpace = inlineTextBuilder.stashedSpace;
        let anyMatch = false;
        let m = wordOrNewlineRe.exec(text);
        if (m) {
          anyMatch = true;
          if (m[0] === '\n') {
            inlineTextBuilder.startNewLine();
          } else if (previouslyStashedSpace || this.testLeadingWhitespace(text)) {
            inlineTextBuilder.pushWord(transform(m[0]));
          } else {
            inlineTextBuilder.concatWord(transform(m[0]));
          }
          while ((m = wordOrNewlineRe.exec(text)) !== null) {
            if (m[0] === '\n') {
              inlineTextBuilder.startNewLine();
            } else {
              inlineTextBuilder.pushWord(transform(m[0]));
            }
          }
        }
        inlineTextBuilder.stashedSpace = (previouslyStashedSpace && !anyMatch) || (this.testTrailingWhitespace(text));
        // No need to stash a space in case last added item was a new line,
        // but that won't affect anything later anyway.
      };

    } else {

      const wordRe = new RegExp(`[^${whitespaceCodes}]+`, 'g');

      this.shrinkWrapAdd = function (text, inlineTextBuilder, transform = (str => str)) {
        if (!text) { return; }
        const previouslyStashedSpace = inlineTextBuilder.stashedSpace;
        let anyMatch = false;
        let m = wordRe.exec(text);
        if (m) {
          anyMatch = true;
          if (previouslyStashedSpace || this.testLeadingWhitespace(text)) {
            inlineTextBuilder.pushWord(transform(m[0]));
          } else {
            inlineTextBuilder.concatWord(transform(m[0]));
          }
          while ((m = wordRe.exec(text)) !== null) {
            inlineTextBuilder.pushWord(transform(m[0]));
          }
        }
        inlineTextBuilder.stashedSpace = (previouslyStashedSpace && !anyMatch) || this.testTrailingWhitespace(text);
      };

    }
  }

  /**
   * Test whether the given text starts with HTML whitespace character.
   *
   * @param   { string }  text  The string to test.
   * @returns { boolean }
   */
  testLeadingWhitespace (text) {
    return this.leadingWhitespaceRe.test(text);
  }

  /**
   * Test whether the given text ends with HTML whitespace character.
   *
   * @param   { string }  text  The string to test.
   * @returns { boolean }
   */
  testTrailingWhitespace (text) {
    return this.trailingWhitespaceRe.test(text);
  }

  /**
   * Test whether the given text contains any non-whitespace characters.
   *
   * @param   { string }  text  The string to test.
   * @returns { boolean }
   */
  testContainsWords (text) {
    return !this.allWhitespaceOrEmptyRe.test(text);
  }

  /**
   * Return the number of newlines if there are no words.
   *
   * If any word is found then return zero regardless of the actual number of newlines.
   *
   * @param   { string }  text  Input string.
   * @returns { number }
   */
  countNewlinesNoWords (text) {
    this.newlineOrNonWhitespaceRe.lastIndex = 0;
    let counter = 0;
    let match;
    while ((match = this.newlineOrNonWhitespaceRe.exec(text)) !== null) {
      if (match[0] === '\n') {
        counter++;
      } else {
        return 0;
      }
    }
    return counter;
  }

};

var whitespaceProcessor = { WhitespaceProcessor: WhitespaceProcessor$1 };

const { trimCharacter: trimCharacter$1 } = helper;
// eslint-disable-next-line no-unused-vars
const { StackItem, BlockStackItem, TableCellStackItem, TableRowStackItem, TableStackItem, TransformerStackItem }
  = stackItem;
const { tableToString } = tablePrinter;
const { WhitespaceProcessor } = whitespaceProcessor;

// eslint-disable-next-line import/no-unassigned-import



/**
 * Helps to build text from inline and block elements.
 *
 * @class BlockTextBuilder
 */
let BlockTextBuilder$1 = class BlockTextBuilder {

  /**
   * Creates an instance of BlockTextBuilder.
   *
   * @param { Options } options HtmlToText options.
   * @param { Picker<DomNode, TagDefinition> } picker Selectors decision tree picker.
   */
  constructor (options, picker) {
    this.options = options;
    this.picker = picker;
    this.whitespaceProcessor = new WhitespaceProcessor(options);
    /** @type { StackItem } */
    this._stackItem = new BlockStackItem(options);
    /** @type { TransformerStackItem } */
    this._wordTransformer = undefined;
  }

  /**
   * Put a word-by-word transform function onto the transformations stack.
   *
   * Mainly used for uppercasing. Can be bypassed to add unformatted text such as URLs.
   *
   * Word transformations applied before wrapping.
   *
   * @param { (str: string) => string } wordTransform Word transformation function.
   */
  pushWordTransform (wordTransform) {
    this._wordTransformer = new TransformerStackItem(this._wordTransformer, wordTransform);
  }

  /**
   * Remove a function from the word transformations stack.
   *
   * @returns { (str: string) => string } A function that was removed.
   */
  popWordTransform () {
    if (!this._wordTransformer) { return undefined; }
    const transform = this._wordTransformer.transform;
    this._wordTransformer = this._wordTransformer.next;
    return transform;
  }

  /** @returns { (str: string) => string } */
  _getCombinedWordTransformer () {
    const applyTransformer = (str, transformer) =>
      ((transformer) ? applyTransformer(transformer.transform(str), transformer.next) : str);
    return (str) => applyTransformer(str, this._wordTransformer);
  }

  _popStackItem () {
    const item = this._stackItem;
    this._stackItem = item.next;
    return item;
  }

  /**
   * Add a line break into currently built block.
   */
  addLineBreak () {
    if (!(
      this._stackItem instanceof BlockStackItem
      || this._stackItem instanceof TableCellStackItem
    )) { return; }
    if (this._stackItem.isPre) {
      this._stackItem.rawText += '\n';
    } else {
      this._stackItem.inlineTextBuilder.startNewLine();
    }
  }

  /**
   * Allow to break line in case directly following text will not fit.
   */
  addWordBreakOpportunity () {
    if (
      this._stackItem instanceof BlockStackItem
      || this._stackItem instanceof TableCellStackItem
    ) {
      this._stackItem.inlineTextBuilder.wordBreakOpportunity = true;
    }
  }

  /**
   * Add a node inline into the currently built block.
   *
   * @param { string } str
   * Text content of a node to add.
   *
   * @param { object | boolean } [ optionsObjectOrNoWordTransform ]
   * Object holding the parameters of the operation.
   *
   * Boolean value is deprecated.
   *
   * @param { boolean } [ optionsObjectOrNoWordTransform.noWordTransform = false ]
   * Ignore word transformers if there are any.
   */
  addInline (str, optionsObjectOrNoWordTransform = {}) {
    if (typeof optionsObjectOrNoWordTransform === 'object') {
      this._addInline(str, optionsObjectOrNoWordTransform);
    } else {
      this._addInline(str, { noWordTransform: optionsObjectOrNoWordTransform });
    }
  }

  _addInline (str, { noWordTransform = false } = {}) {
    if (!(
      this._stackItem instanceof BlockStackItem
      || this._stackItem instanceof TableCellStackItem
    )) { return; }

    if (this._stackItem.isPre) {
      this._stackItem.rawText += str;
      return;
    }

    if (
      str.length === 0 || // empty string
      (
        this._stackItem.stashedLineBreaks && // stashed linebreaks make whitespace irrelevant
        !this.whitespaceProcessor.testContainsWords(str) // no words to add
      )
    ) { return; }

    if (this.options.preserveNewlines) {
      const newlinesNumber = this.whitespaceProcessor.countNewlinesNoWords(str);
      if (newlinesNumber > 0) {
        this._stackItem.inlineTextBuilder.startNewLine(newlinesNumber);
        // keep stashedLineBreaks unchanged
        return;
      }
    }

    if (this._stackItem.stashedLineBreaks) {
      this._stackItem.inlineTextBuilder.startNewLine(this._stackItem.stashedLineBreaks);
    }
    this.whitespaceProcessor.shrinkWrapAdd(
      str,
      this._stackItem.inlineTextBuilder,
      (this._wordTransformer && !noWordTransform) ? this._getCombinedWordTransformer() : undefined
    );
    this._stackItem.stashedLineBreaks = 0; // inline text doesn't introduce line breaks
  }

  /**
   * Start building a new block.
   *
   * @param { object | number } [optionsObjectOrLeadingLineBreaks]
   * Object holding the parameters of the block.
   *
   * Number value is deprecated.
   *
   * @param { number }  [optionsObjectOrLeadingLineBreaks.leadingLineBreaks = 1]
   * This block should have at least this number of line breaks to separate if from any preceding block.
   *
   * @param { number }  [optionsObjectOrLeadingLineBreaks.reservedLineLength = 0]
   * Reserve this number of characters on each line for block markup.
   *
   * @param { boolean } [optionsObjectOrLeadingLineBreaks.isPre = false]
   * Should HTML whitespace be preserved inside this block.
   *
   * @param { number }  [reservedLineLength]
   * Deprecated.
   *
   * @param { boolean } [isPre]
   * Deprecated.
   */
  openBlock (optionsObjectOrLeadingLineBreaks = {}, reservedLineLength = undefined, isPre = undefined) {
    if (typeof optionsObjectOrLeadingLineBreaks === 'object') {
      this._openBlock(optionsObjectOrLeadingLineBreaks);
    } else {
      this._openBlock({
        isPre: isPre,
        leadingLineBreaks: optionsObjectOrLeadingLineBreaks,
        reservedLineLength: reservedLineLength,
      });
    }
  }

  _openBlock ({ leadingLineBreaks = 1, reservedLineLength = 0, isPre = false } = {}) {
    const maxLineLength = Math.max(20, this._stackItem.inlineTextBuilder.maxLineLength - reservedLineLength);
    this._stackItem = new BlockStackItem(
      this.options,
      this._stackItem,
      leadingLineBreaks,
      maxLineLength
    );
    if (isPre) { this._stackItem.isPre = true; }
  }

  /**
   * Finalize currently built block, add it's content to the parent block.
   *
   * @param { object | number }         [optionsObjectOrTrailingLineBreaks]
   * Object holding the parameters of the block.
   *
   * Number value is deprecated.
   *
   * @param { number }                  [optionsObjectOrTrailingLineBreaks.trailingLineBreaks = 1]
   * This block should have at least this number of line breaks to separate it from any following block.
   *
   * @param { (str: string) => string } [optionsObjectOrTrailingLineBreaks.blockTransform = undefined]
   * A function to transform the block text before adding to the parent block.
   * This happens after word wrap and should be used in combination with reserved line length
   * in order to keep line lengths correct.
   * Used for whole block markup.
   *
   * @param { (str: string) => string } [blockTransform]
   * Deprecated.
   */
  closeBlock (optionsObjectOrTrailingLineBreaks = {}, blockTransform = undefined) {
    if (typeof optionsObjectOrTrailingLineBreaks === 'object') {
      this._closeBlock(optionsObjectOrTrailingLineBreaks);
    } else {
      this._closeBlock({
        trailingLineBreaks: optionsObjectOrTrailingLineBreaks,
        blockTransform: blockTransform,
      });
    }
  }

  _closeBlock ({ trailingLineBreaks = 1, blockTransform = undefined } = {}) {
    const block = this._popStackItem();
    const blockText = (blockTransform) ? blockTransform(getText(block)) : getText(block);
    addText(this._stackItem, blockText, block.leadingLineBreaks, Math.max(block.stashedLineBreaks, trailingLineBreaks));
  }

  /**
   * Start building a table.
   */
  openTable () {
    this._stackItem = new TableStackItem(this._stackItem);
  }

  /**
   * Start building a table row.
   */
  openTableRow () {
    if (!(this._stackItem instanceof TableStackItem)) {
      throw new Error('Can\'t add table row to something that is not a table! Check the formatter.');
    }
    this._stackItem = new TableRowStackItem(this._stackItem);
  }

  /**
   * Start building a table cell.
   *
   * @param { object | number } [optionsObjectOrMaxColumnWidth = undefined]
   * Object holding the parameters of the cell.
   *
   * Number value is deprecated.
   *
   * @param { number } [optionsObjectOrMaxColumnWidth.maxColumnWidth = undefined]
   * Wrap cell content to this width. Fall back to global wordwrap value if undefined.
   */
  openTableCell (optionsObjectOrMaxColumnWidth = {}) {
    if (typeof optionsObjectOrMaxColumnWidth === 'object') {
      this._openTableCell(optionsObjectOrMaxColumnWidth);
    } else {
      this._openTableCell({ maxColumnWidth: optionsObjectOrMaxColumnWidth });
    }
  }

  _openTableCell ({ maxColumnWidth = undefined } = {}) {
    if (!(this._stackItem instanceof TableRowStackItem)) {
      throw new Error('Can\'t add table cell to something that is not a table row! Check the formatter.');
    }
    this._stackItem = new TableCellStackItem(this.options, this._stackItem, maxColumnWidth);
  }

  /**
   * Finalize currently built table cell and add it to parent table row's cells.
   *
   * @param { object | number } [optionsObjectOrColspan]
   * Object holding the parameters of the cell.
   *
   * Number value is deprecated.
   *
   * @param { number } [optionsObjectOrColspan.colspan = 1] How many columns this cell should occupy.
   * @param { number } [optionsObjectOrColspan.rowspan = 1] How many rows this cell should occupy.
   *
   * @param { number } [rowspan] Deprecated.
   */
  closeTableCell (optionsObjectOrColspan = {}, rowspan = undefined) {
    if (typeof optionsObjectOrColspan === 'object') {
      this._closeTableCell(optionsObjectOrColspan);
    } else {
      this._closeTableCell({
        colspan: optionsObjectOrColspan,
        rowspan: rowspan,
      });
    }
  }

  _closeTableCell ({ colspan = 1, rowspan = 1 } = {}) {
    const cell = this._popStackItem();
    const text = trimCharacter$1(getText(cell), '\n');
    cell.next.cells.push({ colspan: colspan, rowspan: rowspan, text: text });
  }

  /**
   * Finalize currently built table row and add it to parent table's rows.
   */
  closeTableRow () {
    const row = this._popStackItem();
    row.next.rows.push(row.cells);
  }

  /**
   * Finalize currently built table and add the rendered text to the parent block.
   *
   * @param { object | number } [optionsObjectOrColSpacing]
   * Object holding the parameters of the table.
   *
   * Number value is deprecated.
   *
   * @param { number } [optionsObjectOrColSpacing.colSpacing = 3]
   * Number of spaces between table columns.
   *
   * @param { number } [optionsObjectOrColSpacing.rowSpacing = 0]
   * Number of empty lines between table rows.
   *
   * @param { number } [optionsObjectOrColSpacing.leadingLineBreaks = 2]
   * This table should have at least this number of line breaks to separate if from any preceding block.
   *
   * @param { number } [optionsObjectOrColSpacing.trailingLineBreaks = 2]
   * This table should have at least this number of line breaks to separate it from any following block.
   *
   * @param { number } [rowSpacing]
   * Deprecated.
   *
   * @param { number } [leadingLineBreaks]
   * Deprecated.
   *
   * @param { number } [trailingLineBreaks]
   * Deprecated.
   */
  closeTable (
    optionsObjectOrColSpacing = {},
    rowSpacing = undefined,
    leadingLineBreaks = undefined,
    trailingLineBreaks = undefined
  ) {
    if (typeof optionsObjectOrColSpacing === 'object') {
      this._closeTable(optionsObjectOrColSpacing);
    } else {
      this._closeTable({
        colSpacing: optionsObjectOrColSpacing,
        leadingLineBreaks: leadingLineBreaks,
        rowSpacing: rowSpacing,
        trailingLineBreaks: trailingLineBreaks
      });
    }
  }

  _closeTable ({ colSpacing = 3, rowSpacing = 0, leadingLineBreaks = 2, trailingLineBreaks = 2 } = {}) {
    const table = this._popStackItem();
    const output = tableToString(table.rows, rowSpacing, colSpacing);
    if (output) {
      addText(this._stackItem, output, leadingLineBreaks, trailingLineBreaks);
    }
  }

  /**
   * Return the rendered text content of this builder.
   *
   * @returns { string }
   */
  toString () {
    return getText(this._stackItem.getRoot());
    // There should only be the root item if everything is closed properly.
  }

};

function getText (stackItem) {
  if (!(
    stackItem instanceof BlockStackItem
    || stackItem instanceof TableCellStackItem
  )) {
    throw new Error('Only blocks and table cells can be requested for text contents.');
  }
  return (stackItem.inlineTextBuilder.isEmpty())
    ? stackItem.rawText
    : stackItem.rawText + stackItem.inlineTextBuilder.toString();
}

function addText (stackItem, text, leadingLineBreaks, trailingLineBreaks) {
  if (!(
    stackItem instanceof BlockStackItem
    || stackItem instanceof TableCellStackItem
  )) {
    throw new Error('Only blocks and table cells can contain text.');
  }
  const parentText = getText(stackItem);
  const lineBreaks = Math.max(stackItem.stashedLineBreaks, leadingLineBreaks);
  stackItem.inlineTextBuilder.clear();
  if (parentText) {
    stackItem.rawText = parentText + '\n'.repeat(lineBreaks) + text;
  } else {
    stackItem.rawText = text;
    stackItem.leadingLineBreaks = lineBreaks;
  }
  stackItem.stashedLineBreaks = trailingLineBreaks;
}

var blockTextBuilder = { BlockTextBuilder: BlockTextBuilder$1 };

const he$1 = heExports;

const { get: get$1, numberToLetterSequence, numberToRoman, splitClassesAndIds, trimCharacter } = helper;

// eslint-disable-next-line import/no-unassigned-import



/**
 * Dummy formatter that discards the input and does nothing.
 *
 * @type { FormatCallback }
 */
function formatSkip (elem, walk, builder, formatOptions) {
  /* do nothing */
}

/**
 * Process an inline-level element.
 *
 * @type { FormatCallback }
 */
function formatInline (elem, walk, builder, formatOptions) {
  walk(elem.children, builder);
}

/**
 * Process a block-level container.
 *
 * @type { FormatCallback }
 */
function formatBlock (elem, walk, builder, formatOptions) {
  builder.openBlock({ leadingLineBreaks: formatOptions.leadingLineBreaks });
  walk(elem.children, builder);
  builder.closeBlock({ trailingLineBreaks: formatOptions.trailingLineBreaks });
}

/**
 * Process a line-break.
 *
 * @type { FormatCallback }
 */
function formatLineBreak (elem, walk, builder, formatOptions) {
  builder.addLineBreak();
}

/**
 * Process a `wbk` tag (word break opportunity).
 *
 * @type { FormatCallback }
 */
function formatWbr (elem, walk, builder, formatOptions) {
  builder.addWordBreakOpportunity();
}

/**
 * Process a horizontal line.
 *
 * @type { FormatCallback }
 */
function formatHorizontalLine (elem, walk, builder, formatOptions) {
  builder.openBlock({ leadingLineBreaks: formatOptions.leadingLineBreaks || 2 });
  builder.addInline('-'.repeat(formatOptions.length || builder.options.wordwrap || 40));
  builder.closeBlock({ trailingLineBreaks: formatOptions.trailingLineBreaks || 2 });
}

/**
 * Process a paragraph.
 *
 * @type { FormatCallback }
 */
function formatParagraph (elem, walk, builder, formatOptions) {
  builder.openBlock({ leadingLineBreaks: formatOptions.leadingLineBreaks || 2 });
  walk(elem.children, builder);
  builder.closeBlock({ trailingLineBreaks: formatOptions.trailingLineBreaks || 2 });
}

/**
 * Process a preformatted content.
 *
 * @type { FormatCallback }
 */
function formatPre (elem, walk, builder, formatOptions) {
  builder.openBlock({
    isPre: true,
    leadingLineBreaks: formatOptions.leadingLineBreaks || 2
  });
  walk(elem.children, builder);
  builder.closeBlock({ trailingLineBreaks: formatOptions.trailingLineBreaks || 2 });
}

/**
 * Process a heading.
 *
 * @type { FormatCallback }
 */
function formatHeading (elem, walk, builder, formatOptions) {
  builder.openBlock({ leadingLineBreaks: formatOptions.leadingLineBreaks || 2 });
  if (formatOptions.uppercase !== false) {
    builder.pushWordTransform(str => str.toUpperCase());
    walk(elem.children, builder);
    builder.popWordTransform();
  } else {
    walk(elem.children, builder);
  }
  builder.closeBlock({ trailingLineBreaks: formatOptions.trailingLineBreaks || 2 });
}

/**
 * Process a blockquote.
 *
 * @type { FormatCallback }
 */
function formatBlockquote (elem, walk, builder, formatOptions) {
  builder.openBlock({
    leadingLineBreaks: formatOptions.leadingLineBreaks || 2,
    reservedLineLength: 2
  });
  walk(elem.children, builder);
  builder.closeBlock({
    trailingLineBreaks: formatOptions.trailingLineBreaks || 2,
    blockTransform: str => ((formatOptions.trimEmptyLines !== false) ? trimCharacter(str, '\n') : str)
      .split('\n')
      .map(line => '> ' + line)
      .join('\n')
  });
}

function withBrackets (str, brackets) {
  if (!brackets) { return str; }

  const lbr = (typeof brackets[0] === 'string')
    ? brackets[0]
    : '[';
  const rbr = (typeof brackets[1] === 'string')
    ? brackets[1]
    : ']';
  return lbr + str + rbr;
}

/**
 * Process an image.
 *
 * @type { FormatCallback }
 */
function formatImage (elem, walk, builder, formatOptions) {
  const attribs = elem.attribs || {};
  const alt = (attribs.alt)
    ? he$1.decode(attribs.alt, builder.options.decodeOptions)
    : '';
  const src = (!attribs.src)
    ? ''
    : (formatOptions.baseUrl && attribs.src.indexOf('/') === 0)
      ? formatOptions.baseUrl + attribs.src
      : attribs.src;
  const text = (!src)
    ? alt
    : (!alt)
      ? withBrackets(src, formatOptions.linkBrackets)
      : alt + ' ' + withBrackets(src, formatOptions.linkBrackets);

  builder.addInline(text);
}

/**
 * Process an anchor.
 *
 * @type { FormatCallback }
 */
function formatAnchor (elem, walk, builder, formatOptions) {
  function getHref () {
    if (formatOptions.ignoreHref) { return ''; }
    if (!elem.attribs || !elem.attribs.href) { return ''; }
    let href = elem.attribs.href.replace(/^mailto:/, '');
    if (formatOptions.noAnchorUrl && href[0] === '#') { return ''; }
    href = (formatOptions.baseUrl && href[0] === '/')
      ? formatOptions.baseUrl + href
      : href;
    return he$1.decode(href, builder.options.decodeOptions);
  }
  const href = getHref();
  if (!href) {
    walk(elem.children, builder);
  } else {
    let text = '';
    builder.pushWordTransform(
      str => {
        if (str) { text += str; }
        return str;
      }
    );
    walk(elem.children, builder);
    builder.popWordTransform();

    const hideSameLink = formatOptions.hideLinkHrefIfSameAsText && href === text;
    if (!hideSameLink) {
      builder.addInline(
        (!text)
          ? href
          : ' ' + withBrackets(href, formatOptions.linkBrackets),
        { noWordTransform: true }
      );
    }
  }
}

/**
 * @param { DomNode }           elem               List items with their prefixes.
 * @param { RecursiveCallback } walk               Recursive callback to process child nodes.
 * @param { BlockTextBuilder }  builder            Passed around to accumulate output text.
 * @param { FormatOptions }     formatOptions      Options specific to a formatter.
 * @param { () => string }      nextPrefixCallback Function that returns increasing index each time it is called.
 */
function formatList (elem, walk, builder, formatOptions, nextPrefixCallback) {
  const isNestedList = get$1(elem, ['parent', 'name']) === 'li';

  // With Roman numbers, index length is not as straightforward as with Arabic numbers or letters,
  // so the dumb length comparison is the most robust way to get the correct value.
  let maxPrefixLength = 0;
  const listItems = (elem.children || [])
    // it might be more accurate to check only for html spaces here, but no significant benefit
    .filter(child => child.type !== 'text' || !/^\s*$/.test(child.data))
    .map(function (child) {
      if (child.name !== 'li') {
        return { node: child, prefix: '' };
      }
      const prefix = (isNestedList)
        ? nextPrefixCallback().trimStart()
        : nextPrefixCallback();
      if (prefix.length > maxPrefixLength) { maxPrefixLength = prefix.length; }
      return { node: child, prefix: prefix };
    });
  if (!listItems.length) { return; }

  const reservedLineLength = maxPrefixLength;
  const spacing = '\n' + ' '.repeat(reservedLineLength);
  builder.openBlock({ leadingLineBreaks: isNestedList ? 1 : (formatOptions.leadingLineBreaks || 2) });
  for (const { node, prefix } of listItems) {
    builder.openBlock({
      leadingLineBreaks: 1,
      reservedLineLength: reservedLineLength
    });
    walk([node], builder);
    builder.closeBlock({
      trailingLineBreaks: 1,
      blockTransform: str => prefix + ' '.repeat(reservedLineLength - prefix.length) + str.replace(/\n/g, spacing)
    });
  }
  builder.closeBlock({ trailingLineBreaks: isNestedList ? 1 : (formatOptions.trailingLineBreaks || 2) });
}

/**
 * Process an unordered list.
 *
 * @type { FormatCallback }
 */
function formatUnorderedList (elem, walk, builder, formatOptions) {
  const prefix = formatOptions.itemPrefix || ' * ';
  return formatList(elem, walk, builder, formatOptions, () => prefix);
}

/**
 * Process an ordered list.
 *
 * @type { FormatCallback }
 */
function formatOrderedList (elem, walk, builder, formatOptions) {
  let nextIndex = Number(elem.attribs.start || '1');
  const indexFunction = getOrderedListIndexFunction(elem.attribs.type);
  const nextPrefixCallback = () => ' ' + indexFunction(nextIndex++) + '. ';
  return formatList(elem, walk, builder, formatOptions, nextPrefixCallback);
}

/**
 * Return a function that can be used to generate index markers of a specified format.
 *
 * @param   { string } [olType='1'] Marker type.
 * @returns { (i: number) => string }
 */
function getOrderedListIndexFunction (olType = '1') {
  switch (olType) {
    case 'a': return (i) => numberToLetterSequence(i, 'a');
    case 'A': return (i) => numberToLetterSequence(i, 'A');
    case 'i': return (i) => numberToRoman(i).toLowerCase();
    case 'I': return (i) => numberToRoman(i);
    case '1':
    default: return (i) => (i).toString();
  }
}

function isDataTable (attr, tables) {
  if (tables === true) { return true; }
  if (!attr) { return false; }

  const { classes, ids } = splitClassesAndIds(tables);
  const attrClasses = (attr['class'] || '').split(' ');
  const attrIds = (attr['id'] || '').split(' ');

  return attrClasses.some(x => classes.includes(x)) || attrIds.some(x => ids.includes(x));
}

/**
 * Process a table (either as a container or as a data table, depending on options).
 *
 * @type { FormatCallback }
 */
function formatTable (elem, walk, builder, formatOptions) {
  return isDataTable(elem.attribs, builder.options.tables)
    ? formatDataTable(elem, walk, builder, formatOptions)
    : formatBlock(elem, walk, builder, formatOptions);
}

/**
 * Process a data table.
 *
 * @type { FormatCallback }
 */
function formatDataTable (elem, walk, builder, formatOptions) {
  builder.openTable();
  elem.children.forEach(walkTable);
  builder.closeTable({
    colSpacing: formatOptions.colSpacing,
    leadingLineBreaks: formatOptions.leadingLineBreaks,
    rowSpacing: formatOptions.rowSpacing,
    trailingLineBreaks: formatOptions.trailingLineBreaks
  });

  function formatCell (cellNode) {
    const colspan = +get$1(cellNode, ['attribs', 'colspan']) || 1;
    const rowspan = +get$1(cellNode, ['attribs', 'rowspan']) || 1;
    builder.openTableCell({ maxColumnWidth: formatOptions.maxColumnWidth });
    walk(cellNode.children, builder);
    builder.closeTableCell({ colspan: colspan, rowspan: rowspan });
  }

  function walkTable (elem) {
    if (elem.type !== 'tag') { return; }

    const formatHeaderCell = (formatOptions.uppercaseHeaderCells !== false)
      ? (cellNode) => {
        builder.pushWordTransform(str => str.toUpperCase());
        formatCell(cellNode);
        builder.popWordTransform();
      }
      : formatCell;

    switch (elem.name) {
      case 'thead':
      case 'tbody':
      case 'tfoot':
      case 'center':
        elem.children.forEach(walkTable);
        return;

      case 'tr': {
        builder.openTableRow();
        for (const childOfTr of elem.children) {
          if (childOfTr.type !== 'tag') { continue; }
          switch (childOfTr.name) {
            case 'th': {
              formatHeaderCell(childOfTr);
              break;
            }
            case 'td': {
              formatCell(childOfTr);
              break;
            }
              // do nothing
          }
        }
        builder.closeTableRow();
        break;
      }
        // do nothing
    }
  }
}

var formatter = {
  anchor: formatAnchor,
  block: formatBlock,
  blockquote: formatBlockquote,
  dataTable: formatDataTable,
  heading: formatHeading,
  horizontalLine: formatHorizontalLine,
  image: formatImage,
  inline: formatInline,
  lineBreak: formatLineBreak,
  orderedList: formatOrderedList,
  paragraph: formatParagraph,
  pre: formatPre,
  skip: formatSkip,
  table: formatTable,
  unorderedList: formatUnorderedList,
  wbr: formatWbr
};

const { hp2Builder } = hp2Builder$2;
const merge = cjs;
const he = heExports;
const htmlparser = lib$3;
const selderee = selderee$2;

const { BlockTextBuilder } = blockTextBuilder;
const defaultFormatters = formatter;
const { limitedDepthRecursive, mergeDuplicatesPreferLast, get } = helper;

// eslint-disable-next-line import/no-unassigned-import



/**
 * Default options.
 *
 * @constant
 * @type { Options }
 * @default
 * @private
 */
const DEFAULT_OPTIONS = {
  baseElements: {
    selectors: [ 'body' ],
    orderBy: 'selectors', // 'selectors' | 'occurrence'
    returnDomByDefault: true
  },
  decodeOptions: {
    isAttributeValue: false,
    strict: false
  },
  formatters: {},
  limits: {
    ellipsis: '...',
    maxBaseElements: undefined,
    maxChildNodes: undefined,
    maxDepth: undefined,
    maxInputLength: (1 << 24) // 16_777_216
  },
  longWordSplit: {
    forceWrapOnLimit: false,
    wrapCharacters: []
  },
  preserveNewlines: false,
  selectors: [
    { selector: '*', format: 'inline' },
    {
      selector: 'a',
      format: 'anchor',
      options: {
        baseUrl: null,
        hideLinkHrefIfSameAsText: false,
        ignoreHref: false,
        linkBrackets: ['[', ']'],
        noAnchorUrl: true
      }
    },
    { selector: 'article', format: 'block' },
    { selector: 'aside', format: 'block' },
    {
      selector: 'blockquote',
      format: 'blockquote',
      options: { leadingLineBreaks: 2, trailingLineBreaks: 2, trimEmptyLines: true }
    },
    { selector: 'br', format: 'lineBreak' },
    { selector: 'div', format: 'block' },
    { selector: 'footer', format: 'block' },
    { selector: 'form', format: 'block' },
    { selector: 'h1', format: 'heading', options: { leadingLineBreaks: 3, trailingLineBreaks: 2, uppercase: true } },
    { selector: 'h2', format: 'heading', options: { leadingLineBreaks: 3, trailingLineBreaks: 2, uppercase: true } },
    { selector: 'h3', format: 'heading', options: { leadingLineBreaks: 3, trailingLineBreaks: 2, uppercase: true } },
    { selector: 'h4', format: 'heading', options: { leadingLineBreaks: 2, trailingLineBreaks: 2, uppercase: true } },
    { selector: 'h5', format: 'heading', options: { leadingLineBreaks: 2, trailingLineBreaks: 2, uppercase: true } },
    { selector: 'h6', format: 'heading', options: { leadingLineBreaks: 2, trailingLineBreaks: 2, uppercase: true } },
    { selector: 'header', format: 'block' },
    {
      selector: 'hr',
      format: 'horizontalLine',
      options: { leadingLineBreaks: 2, length: undefined, trailingLineBreaks: 2 }
    },
    {
      selector: 'img',
      format: 'image',
      options: { baseUrl: null, linkBrackets: ['[', ']'] }
    },
    { selector: 'main', format: 'block' },
    { selector: 'nav', format: 'block' },
    {
      selector: 'ol',
      format: 'orderedList',
      options: { leadingLineBreaks: 2, trailingLineBreaks: 2 }
    },
    { selector: 'p', format: 'paragraph', options: { leadingLineBreaks: 2, trailingLineBreaks: 2 } },
    { selector: 'pre', format: 'pre', options: { leadingLineBreaks: 2, trailingLineBreaks: 2 } },
    { selector: 'section', format: 'block' },
    {
      selector: 'table',
      format: 'table',
      options: {
        colSpacing: 3,
        leadingLineBreaks: 2,
        maxColumnWidth: 60,
        rowSpacing: 0,
        trailingLineBreaks: 2,
        uppercaseHeaderCells: true
      }
    },
    {
      selector: 'ul',
      format: 'unorderedList',
      options: { itemPrefix: ' * ', leadingLineBreaks: 2, trailingLineBreaks: 2 }
    },
    { selector: 'wbr', format: 'wbr' },
  ],
  tables: [], // deprecated
  whitespaceCharacters: ' \t\r\n\f\u200b',
  wordwrap: 80
};

const concatMerge = (acc, src, options) => [...acc, ...src];
const overwriteMerge = (acc, src, options) => [...src];
const selectorsMerge = (acc, src, options) => (
  (acc.some(s => typeof s === 'object'))
    ? concatMerge(acc, src) // selectors
    : overwriteMerge(acc, src) // baseElements.selectors
);

/**
 * Preprocess options, compile selectors into a decision tree,
 * return a function intended for batch processing.
 *
 * @param   { Options } [options = {}]   HtmlToText options.
 * @returns { (html: string) => string } Pre-configured converter function.
 * @static
 */
function compile (options = {}) {
  options = merge(
    DEFAULT_OPTIONS,
    options,
    {
      arrayMerge: overwriteMerge,
      customMerge: (key) => ((key === 'selectors') ? selectorsMerge : undefined)
    }
  );
  options.formatters = Object.assign({}, defaultFormatters, options.formatters);
  options.selectors = mergeDuplicatesPreferLast(options.selectors, (s => s.selector));

  handleDeprecatedOptions(options);

  const selectorsWithoutFormat = options.selectors.filter(s => !s.format);
  if (selectorsWithoutFormat.length) {
    throw new Error(
      'Following selectors have no specified format: ' +
      selectorsWithoutFormat.map(s => `\`${s.selector}\``).join(', ')
    );
  }
  const picker = new selderee.DecisionTree(
    options.selectors.map(s => [s.selector, s])
  ).build(hp2Builder);

  const baseSelectorsPicker = new selderee.DecisionTree(
    options.baseElements.selectors.map((s, i) => [s, i + 1])
  ).build(hp2Builder);
  function findBaseElements (dom) {
    return findBases(dom, options, baseSelectorsPicker);
  }

  const limitedWalk = limitedDepthRecursive(
    options.limits.maxDepth,
    recursiveWalk,
    function (dom, builder) {
      builder.addInline(options.limits.ellipsis || '');
    }
  );

  return function (html) {
    return process(html, options, picker, findBaseElements, limitedWalk);
  };
}

/**
 * Convert given HTML according to preprocessed options.
 *
 * @param { string } html HTML content to convert.
 * @param { Options } options HtmlToText options (preprocessed).
 * @param { Picker<DomNode, TagDefinition> } picker
 * Tag definition picker for DOM nodes processing.
 * @param { (dom: DomNode[]) => DomNode[] } findBaseElements
 * Function to extract elements from HTML DOM
 * that will only be present in the output text.
 * @param { RecursiveCallback } walk Recursive callback.
 * @returns { string }
 */
function process (html, options, picker, findBaseElements, walk) {
  const maxInputLength = options.limits.maxInputLength;
  if (maxInputLength && html && html.length > maxInputLength) {
    console.warn(
      `Input length ${html.length} is above allowed limit of ${maxInputLength}. Truncating without ellipsis.`
    );
    html = html.substring(0, maxInputLength);
  }

  const handler = new htmlparser.DomHandler();
  new htmlparser.Parser(handler, { decodeEntities: false }).parseComplete(html);

  const bases = findBaseElements(handler.dom);
  const builder = new BlockTextBuilder(options, picker);
  walk(bases, builder);
  return builder.toString();
}

/**
 * Convert given HTML content to plain text string.
 *
 * @param   { string }  html           HTML content to convert.
 * @param   { Options } [options = {}] HtmlToText options.
 * @returns { string }                 Plain text string.
 * @static
 *
 * @example
 * const { convert } = require('html-to-text');
 * const text = convert('<h1>Hello World</h1>', {
 *   wordwrap: 130
 * });
 * console.log(text); // HELLO WORLD
 */
function convert (html, options = {}) {
  return compile(options)(html);
}

/**
 * Map previously existing and now deprecated options to the new options layout.
 * This is a subject for cleanup in major releases.
 *
 * @param { Options } options HtmlToText options.
 */
function handleDeprecatedOptions (options) {
  const selectorDefinitions = options.selectors;

  if (options.tags) {
    const tagDefinitions = Object.entries(options.tags).map(
      ([selector, definition]) => ({ ...definition, selector: selector || '*' })
    );
    selectorDefinitions.push(...tagDefinitions);
  }

  function set (obj, path, value) {
    const valueKey = path.pop();
    for (const key of path) {
      let nested = obj[key];
      if (!nested) {
        nested = {};
        obj[key] = nested;
      }
      obj = nested;
    }
    obj[valueKey] = value;
  }

  function copyFormatterOption (source, format, target) {
    if (options[source] === undefined) { return; }
    for (const definition of selectorDefinitions) {
      if (definition.format === format) {
        set(definition, ['options', target], options[source]);
      }
    }
  }

  copyFormatterOption('hideLinkHrefIfSameAsText', 'anchor', 'hideLinkHrefIfSameAsText');
  copyFormatterOption('ignoreHref', 'anchor', 'ignoreHref');
  copyFormatterOption('linkHrefBaseUrl', 'anchor', 'baseUrl');
  copyFormatterOption('noAnchorUrl', 'anchor', 'noAnchorUrl');
  copyFormatterOption('noLinkBrackets', 'anchor', 'noLinkBrackets');

  copyFormatterOption('linkHrefBaseUrl', 'image', 'baseUrl');

  copyFormatterOption('unorderedListItemPrefix', 'unorderedList', 'itemPrefix');

  copyFormatterOption('uppercaseHeadings', 'heading', 'uppercase');
  copyFormatterOption('uppercaseHeadings', 'table', 'uppercaseHeadings');
  copyFormatterOption('uppercaseHeadings', 'dataTable', 'uppercaseHeadings');

  if (options['ignoreImage']) {
    for (const definition of selectorDefinitions) {
      if (definition.format === 'image') {
        definition.format = 'skip';
      }
    }
  }

  if (options['singleNewLineParagraphs']) {
    for (const definition of selectorDefinitions) {
      if (definition.format === 'paragraph' || definition.format === 'pre') {
        set(definition, ['options', 'leadingLineBreaks'], 1);
        set(definition, ['options', 'trailingLineBreaks'], 1);
      }
    }
  }

  if (options['baseElement']) {
    const baseElement = options['baseElement'];
    set(
      options,
      ['baseElements', 'selectors'],
      (Array.isArray(baseElement) ? baseElement : [baseElement])
    );
  }
  if (options['returnDomByDefault'] !== undefined) {
    set(options, ['baseElements', 'returnDomByDefault'], options['returnDomByDefault']);
  }

  for (const definition of selectorDefinitions) {
    if (definition.format === 'anchor' && get(definition, ['options', 'noLinkBrackets'])) {
      set(definition, ['options', 'linkBrackets'], false);
    }
  }
}

function findBases (dom, options, baseSelectorsPicker) {
  const results = [];

  function recursiveWalk (walk, /** @type { DomNode[] } */ dom) {
    dom = dom.slice(0, options.limits.maxChildNodes);
    for (const elem of dom) {
      if (elem.type !== 'tag') {
        continue;
      }
      const pickedSelectorIndex = baseSelectorsPicker.pick1(elem);
      if (pickedSelectorIndex > 0) {
        results.push({ selectorIndex: pickedSelectorIndex, element: elem });
      } else if (elem.children) {
        walk(elem.children);
      }
      if (results.length >= options.limits.maxBaseElements) {
        return;
      }
    }
  }

  const limitedWalk = limitedDepthRecursive(
    options.limits.maxDepth,
    recursiveWalk
  );
  limitedWalk(dom);

  if (options.baseElements.orderBy !== 'occurrence') { // 'selectors'
    results.sort((a, b) => a.selectorIndex - b.selectorIndex);
  }
  return (options.baseElements.returnDomByDefault && results.length === 0)
    ? dom
    : results.map(x => x.element);
}

/**
 * Function to walk through DOM nodes and accumulate their string representations.
 *
 * @param   { RecursiveCallback } walk    Recursive callback.
 * @param   { DomNode[] }         [dom]   Nodes array to process.
 * @param   { BlockTextBuilder }  builder Passed around to accumulate output text.
 * @private
 */
function recursiveWalk (walk, dom, builder) {
  if (!dom) { return; }

  const options = builder.options;

  const tooManyChildNodes = dom.length > options.limits.maxChildNodes;
  if (tooManyChildNodes) {
    dom = dom.slice(0, options.limits.maxChildNodes);
    dom.push({
      data: options.limits.ellipsis,
      type: 'text'
    });
  }

  for (const elem of dom) {
    switch (elem.type) {
      case 'text': {
        builder.addInline(he.decode(elem.data, options.decodeOptions));
        break;
      }
      case 'tag': {
        const tagDefinition = builder.picker.pick1(elem);
        const format = options.formatters[tagDefinition.format];
        format(elem, walk, builder, tagDefinition.options || {});
        break;
      }
    }
  }

  return;
}

/**
 * @deprecated Use `{ convert }` function instead!
 * @see convert
 *
 * @param   { string }  html           HTML content to convert.
 * @param   { Options } [options = {}] HtmlToText options.
 * @returns { string }                 Plain text string.
 * @static
 */
const fromString = (html, options = {}) => convert(html, options);

var htmlToText$1 = {
  compile: compile,
  convert: convert,
  fromString: fromString,
  htmlToText: convert
};

var htmlToText = htmlToText$1;

const htmlToString = (content) => {
    try {
        return htmlToText.htmlToText(content, {
            tables: ["*"],
            selectors: [{ selector: "span.vc-comment-vars", format: "skip" }],
        });
    }
    catch (e) {
        console.error(e);
        return null;
    }
};

class IndexDataProvider {
    constructor(_pCatalogs, _storage) {
        this._pCatalogs = _pCatalogs;
        this._storage = _storage;
        this._pCatalogs.onUpdate(this._getOnChangeRule());
    }
    getCatalogValue(catalogName) {
        return __awaiter(this, void 0, void 0, function* () {
            if (yield this._storage.exists(catalogName))
                return this._getIndexDataFromStorage(catalogName);
            return this._getAndCreateIndexData(catalogName);
        });
    }
    deleteCatalogs() {
        return __awaiter(this, void 0, void 0, function* () {
            const catalogNames = (yield this._pCatalogs.getAll()).map((c) => c.getName());
            return Promise.all(catalogNames.map((catalogName) => __awaiter(this, void 0, void 0, function* () {
                if (yield this._storage.exists(catalogName))
                    yield this._storage.delete(catalogName);
            })));
        });
    }
    _getAndCreateIndexData(catalogName) {
        return __awaiter(this, void 0, void 0, function* () {
            const data = yield this._getIndexData(catalogName);
            yield this._setIndexDataInStorage(catalogName, data);
            return data;
        });
    }
    _getOnChangeRule() {
        return (changeItems) => __awaiter(this, void 0, void 0, function* () {
            const catalogNames = [...new Set(changeItems.map((item) => item.catalog.getName()))];
            yield Promise.all(catalogNames.map((catalogName) => this._getAndCreateIndexData(catalogName)));
        });
    }
    _getIndexData(catalogName) {
        return __awaiter(this, void 0, void 0, function* () {
            const catalog = yield this._pCatalogs.get(catalogName);
            const contentItems = catalog.getArticles();
            const indexDataPromises = contentItems.map((article) => __awaiter(this, void 0, void 0, function* () {
                var _a;
                const content = htmlToString(yield article.getHtmlContent());
                return {
                    path: article.id,
                    pathname: yield article.getPathname(),
                    title: article.getProp("title"),
                    content: content !== null && content !== void 0 ? content : "",
                    tags: (_a = article.getProp("tags")) === null || _a === void 0 ? void 0 : _a.join(" "),
                };
            }));
            return Promise.all(indexDataPromises);
        });
    }
    _setIndexDataInStorage(catalogName, indexData) {
        return this._storage.set(catalogName, JSON.stringify(indexData));
    }
    _getIndexDataFromStorage(catalogName) {
        return __awaiter(this, void 0, void 0, function* () {
            return this._storage.get(catalogName).then((data) => JSON.parse(data));
        });
    }
}

class SearchPlugin extends Plugin {
    get name() {
        return "search";
    }
    onLoad() {
        this.addCommand({
            name: "searchCommand",
            do(_a) {
                return __awaiter(this, arguments, void 0, function* ({ query, catalogName }) {
                    const indexDataProvider = new IndexDataProvider(this.app.catalogs, this.app.storage);
                    const searcher = new FuseSearcher(indexDataProvider);
                    const getCatalogArticleIds = (catalogName) => __awaiter(this, void 0, void 0, function* () {
                        const catalog = yield this.app.catalogs.get(catalogName);
                        return catalog === null || catalog === void 0 ? void 0 : catalog.getArticles().map((a) => a.id);
                    });
                    const getSearchData = (query, catalogName) => __awaiter(this, void 0, void 0, function* () {
                        var _b;
                        const search = (query) => __awaiter(this, void 0, void 0, function* () {
                            const result = yield searcher.search(query, catalogName, yield getCatalogArticleIds(catalogName));
                            return result.length > 0 ? result : null;
                        });
                        return (_b = (yield multiLayoutSearcher(search)(query))) !== null && _b !== void 0 ? _b : [];
                    });
                    const getAllSearchData = (query) => __awaiter(this, void 0, void 0, function* () {
                        var _c;
                        const catalogs = (yield this.app.catalogs.getAll()).map((l) => l.getName());
                        const catalogArticleIds = {};
                        yield Promise.all(catalogs.map((c) => __awaiter(this, void 0, void 0, function* () { return (catalogArticleIds[c] = yield getCatalogArticleIds(c)); })));
                        const search = (query) => __awaiter(this, void 0, void 0, function* () {
                            const result = yield searcher.searchAll(query, catalogArticleIds);
                            return result.length > 0 ? result : null;
                        });
                        return (_c = (yield multiLayoutSearcher(search)(query))) !== null && _c !== void 0 ? _c : [];
                    });
                    return catalogName ? yield getSearchData(query, catalogName) : yield getAllSearchData(query);
                });
            },
        });
        this.addCommand({
            name: "resetSearchData",
            do() {
                return __awaiter(this, void 0, void 0, function* () {
                    const indexDataProvider = new IndexDataProvider(this.app.catalogs, this.app.storage);
                    const searcher = new FuseSearcher(indexDataProvider);
                    yield searcher.resetAllCatalogs();
                });
            },
        });
    }
    onUnload() { }
}

export { SearchPlugin as default };
