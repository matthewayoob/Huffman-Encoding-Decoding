/* File: Huffman.cpp
 *
 * Name: Matthew Ayoob
 *
 * This file details the methods of Huffman Coding and rigorous testing.
 */
#include "Huffman.h"
#include "GUI/SimpleTest.h"
#include "hashmap.h"
#include "vector.h"
#include "priorityqueue.h"
#include <string>
#include "error.h"
#include "map.h"
using namespace std;

//protoyping
Map<char, Vector<Bit>> createBitmap(EncodingTreeNode* tree);

/**
* Deallocates all nodes in a Huffman tree. We've provided this helper function
* to you since we also use it in our test driver and figured you might want
* to make use of it.
*/
void deleteTree(EncodingTreeNode* tree) {
   if (tree != nullptr) {
       deleteTree(tree->zero);
       deleteTree(tree->one);
       delete tree;
   }
}

/**
* Constructs a Huffman coding tree for the given string, using the algorithm
* described in lecture.
*
* If the input string does not contain at least two different characters,
* this function should report an error.
*
* When assembling larger trees out of smaller ones, make sure to set the first
* tree dequeued from the queue to be the zero subtree of the new tree and the
* second tree as the one subtree.
*/
EncodingTreeNode* huffmanTreeFor(const string& str) {
   PriorityQueue<EncodingTreeNode*> priorityQueue;
   Map<char, int> frequencyMap;
   for (char character: str) {
       if (!frequencyMap.containsKey(character)) {
           frequencyMap[character] = 0;
       }
       frequencyMap[character]++;
   }
   if (frequencyMap.size() < 2) {
       error("String needs to contain at least two distinct characters :( ");
   }
   for (char character: frequencyMap) {
       EncodingTreeNode* node = new EncodingTreeNode;
       node->ch = character;
       node->zero = nullptr;
       node->one = nullptr;
       priorityQueue.enqueue(node, frequencyMap[character]);
   }
   while (priorityQueue.size() > 1) {
       int nodeApriority = priorityQueue.peekPriority();
       EncodingTreeNode* nodeA = priorityQueue.dequeue();
       int nodeBpriority = priorityQueue.peekPriority();
       EncodingTreeNode* nodeB = priorityQueue.dequeue();
       EncodingTreeNode* tree = new EncodingTreeNode;
       tree->zero = nodeA;
       tree->one = nodeB;
       priorityQueue.enqueue(tree, (nodeApriority + nodeBpriority));
   }
   return priorityQueue.peek();
}


/**
* Given a Queue<Bit> containing a compressed message and a tree that was used
* to encode those bits, decodes the bits back to the original message.
*
* You can assume the input tree is not null and is not a single character;
* those are edge cases you don't need to handle.
*
* You can assume the input Queue<Bit> is well-formed in that every character
* was encoded correctly, there are no stray bits in the Queue, etc.
*/
string decodeText(Queue<Bit>& bits, EncodingTreeNode* tree) {
    string decodedText = "";
    char character;
    while (!bits.isEmpty()) {
        EncodingTreeNode* node = tree;
        while (node->zero != nullptr && node->one != nullptr && !bits.isEmpty()) {
            Bit bit = bits.dequeue();
            if (bit == 0) {
                node = node->zero;
                character = node->ch;
            }
            else {
                node = node->one;
                character = node->ch;
            }
        }
        decodedText += character;
    }
    return decodedText;
}

/**
* Given a string and a Huffman encoding tree, encodes that text using the tree
* and outputs a Queue<Bit> corresponding to the encoded representation.
*
* The input tree will not be null and will not consist of a single node; these
* are edge cases you don't have to handle. The input tree will contain all
* characters that make up the input string.
*/
Queue<Bit> encodeText(const string& str, EncodingTreeNode* tree) {
   Queue<Bit> encodedText;
   Map<char, Vector<Bit>> bitMap = createBitmap(tree);
   for (char character:str) {
       for(Bit bit: bitMap[character]) {
           encodedText.enqueue(bit);
       }
   }
   return encodedText;
}

/**
* Given a map of char and vector of bytes, a vector of bytes, a vector of bytes,
* and an encoding tree node pointer, this function outputs a map mapping characters
* to their corresponding vector of bytes to be referenced by the encodeText function.
*/
Map<char, Vector<Bit>> createBitmapRec(Map<char, Vector<Bit>>& bitMap, Vector<Bit>& path, EncodingTreeNode* tree) {
   EncodingTreeNode* node = tree;
   // reached end of a path
   if (node->zero == nullptr && node->one == nullptr) {
       bitMap[node->ch] = path;
   }
   if (node->zero != nullptr) {
       path.add(0);
       createBitmapRec(bitMap, path, node->zero);
       path.remove(path.size() - 1);
   }
   if (node->one != nullptr) {
       path.add(1);
       createBitmapRec(bitMap, path, node->one);
       path.remove(path.size() - 1);
   }
   return bitMap;
}

/**
* This wrapper function initializes a map of char to corresponding to
* vectors of bytes and a vector of bytes and a recursive call.
*
* It returns the output of the createBitmapRec, a Map of chars and vectors
* of bytes.
*/
Map<char, Vector<Bit>> createBitmap(EncodingTreeNode* tree) {
   Map<char, Vector<Bit>> bitMap;
   Vector<Bit> path;
   return createBitmapRec(bitMap, path, tree);
}



/**
* Decodes the given Queue<Bit> and Queue<char> into a Huffman coding tree.
*
* You can assume that the input Queues are structured properly in that they
* represent a legal encoding of a tree, that there aren't stray characters
* or bits in them, etc.
*/
EncodingTreeNode* decodeTree(Queue<Bit>& bits, Queue<char>& leaves) {
    EncodingTreeNode* node = new EncodingTreeNode;
    if (bits.size() == 0){
        return node;
    }
    else {
        if (bits.dequeue() == 1) {
            node->zero = decodeTree(bits, leaves);
            node->one= decodeTree(bits, leaves);
        }
        else {
            node->ch = leaves.dequeue();
            node->zero = nullptr;
            node->one = nullptr;
        }
    }
    return node;

}



/**
* Encodes the given Huffman tree as a Queue<Bit> and Queue<char> in the manner
* specified in the assignment handout.
*
* You can assume the input Queues are empty on entry to this function.
*
* You can assume that the Huffman tree provided to you is properly structured,
* in that each internal node will have two children and only the characters in
* the leaves matter, etc.
*/
void encodeTree(EncodingTreeNode* tree, Queue<Bit>& bits, Queue<char>& leaves) {
    if (tree->zero == nullptr) {
        bits.enqueue(0);
        leaves.enqueue(tree->ch);
    }
    else {
        bits.enqueue(1);
        encodeTree(tree->zero, bits, leaves);
        encodeTree(tree->one, bits, leaves);
    }

}

/**
* Compresses the given text string using Huffman coding, producing as output
* a HuffmanResult containing the encoded tree and message.
*
* Your implementation of this function should report an error if there are
* fewer than two distinct characters in the input string.
*/
HuffmanResult compress(const string& text) {
    if (text.size() < 2) {
        error("Too few characters!");
    }
    for (int i = 1; i < text.size(); i++) {
        if (text[i] != text[i -1]) break;
        if (i == text.size() -1) {
            error("Need at least two distinct characters :(");
        }
    }
    EncodingTreeNode* tree = huffmanTreeFor(text);
    Queue<Bit> bitQueue;
    Queue<char> charQueue;
    HuffmanResult res;
    encodeTree(tree, bitQueue, charQueue);
    res.treeBits = bitQueue;
    res.messageBits = encodeText(text, tree);
    res.treeLeaves = charQueue;
    deleteTree(tree);
    return res;
}

/**
* Decompresses the given HuffmanResult and returns the string it represents.
*
* Your implementation may change the file parameter however it sees fit. There
* are no requirements about what it should look like after this function
* returns.
*
* You can assume the input file is well-formed and was created by a correct
* implementation of compress.
*/
string decompress(HuffmanResult& file) {
    if (file.messageBits.size() < 2) {
        error("String needs to contain at least two distinct characters :( ");
    }
    EncodingTreeNode* tree = decodeTree(file.treeBits, file.treeLeaves);
    string res = decodeText(file.messageBits, tree);
    deleteTree(tree);
    return res;
}

/** Exhaustive Testing */

TEST("huffmanTreeFor handles repeated characters and the empty string") {
    EXPECT_ERROR(huffmanTreeFor("AAAAAAAAAAAAAAAA"));
    EXPECT_ERROR(huffmanTreeFor("................"));
    EXPECT_ERROR(huffmanTreeFor("")); // No characters
}


TEST("Handles error for decompress") {
    HuffmanResult file = {
        { 0 },
        {'u'},
        { 0 }
    };

    EXPECT_ERROR(decompress(file));
}


TEST("encodeTree fully handles left path first") {

        EncodingTreeNode* tree = huffmanTreeFor("AB");

        Queue<Bit>  treeBits;
        Queue<char> treeLeaves;

        Queue<Bit>  expectedBits = { 1, 0, 0};
        Queue<char> expectedLeaves = { 'B', 'A' };

        encodeTree(tree, treeBits, treeLeaves);
        EXPECT_EQUAL(treeBits, expectedBits);
        EXPECT_EQUAL(treeLeaves, expectedLeaves);

        deleteTree(tree);
}

TEST("decodeTree properly places child") {
    /* This tree:
     *                 *
     *                / \
     *               O   *
     *                  / \
     *                 *   N

     */
    EncodingTreeNode* tree = new EncodingTreeNode {
            '*',
            new EncodingTreeNode { 'O', nullptr, nullptr },
            new EncodingTreeNode { '*',
            new EncodingTreeNode{ '*', nullptr, nullptr},
            new EncodingTreeNode{ 'N', nullptr, nullptr }
}
};
    Queue<Bit> bits = { 1, 0, 1, 0, 0 };
    EXPECT_NOT_EQUAL(decodeText(bits, tree), "ON");

    deleteTree(tree);
}


TEST("encodeText accounts for improper trees") {
    /* This tree:
     *                 *
     *                / \
     *               *   *
     *
     */
    EncodingTreeNode* tree = new EncodingTreeNode {
           '*',
           new EncodingTreeNode { '*', nullptr, nullptr },
               new EncodingTreeNode { '*', nullptr, nullptr }
            };
            Queue<Bit> tester = { 0, 0, 0 };
        EXPECT_NOT_EQUAL(encodeText("", tree), tester);
        deleteTree(tree);
}

TEST("decodeText handles a misrepresented child") {
    /* All possible characters. */
    string allChars = "";

    /* Try decoding all pairs of adjacent characters. */
    for (size_t i = 0; i < allChars.size(); i += 2) {
        string expected;
        expected += allChars[i];
        expected += allChars[i + 1];
        expected += allChars[i + 1];

        EncodingTreeNode* tree = new EncodingTreeNode {
            ' ',
            new EncodingTreeNode {allChars[i], nullptr, nullptr},
            new EncodingTreeNode {allChars[i + 1], nullptr, nullptr}
        };

        Queue<Bit> bits = { 0, 1, 0 };
        EXPECT_NOT_EQUAL(decodeText(bits, tree), expected);

        deleteTree(tree);
    }
}

#include <limits>

/* Utility function to test if a purported Huffman tree is indeed a Huffman tree.
 * Specifically, this checks that each internal node has either zero or two
 * children. There are other ways you could produce an invalid Huffman tree - for
 * example, by having uninitialized pointers or by linking in a cycle - but we
 * don't test for that here.
 */
bool isEncodingTree(EncodingTreeNode* tree) {
    /* The empty tree is not a Huffman tree. */
    if (tree == nullptr) return false;

    /* If we have one missing child, we should have two missing children. */
    if ((tree->zero == nullptr) != (tree->one == nullptr)) return false;

    /* If we have children at all, they need to be Huffman trees. */
    return tree->zero == nullptr || (isEncodingTree(tree->zero) && isEncodingTree(tree->one));
}

/* Utility function to test if two trees are equal. This is adapted from Section
 * Handout 8 and particularized to Huffman trees.
 */
bool areEqual(EncodingTreeNode* lhs, EncodingTreeNode* rhs) {
    /* Base case: If either is a leaf, both should be. */
    bool lhsLeaf = lhs->zero == nullptr && lhs->one == nullptr;
    bool rhsLeaf = rhs->zero == nullptr && rhs->one == nullptr;
    if (lhsLeaf || rhsLeaf) {
        return lhs->ch == rhs->ch && lhsLeaf == rhsLeaf;
    }

    /* Otherwise, they're both internal nodes. Check that they match. */
    return areEqual(lhs->zero, rhs->zero) && areEqual(lhs->one, rhs->one);
}

/* Utility function to return a string of all possible characters. */
string pangrammaticString() {
    string result;

    char ch = numeric_limits<char>::min();
    result += ch;
    do {
        ch++;
        result += ch;
    } while (ch != numeric_limits<char>::max());

    return result;
}

/* Utility function that makes an inefficient but still valid encoding tree
 * for the given characters.
 */
EncodingTreeNode* strandTreeFor(const string& text, size_t index = 0) {
    if (index == text.size()) error("No characters provided to strandTreeFor.");

    /* We always get a leaf node. */
    EncodingTreeNode* leaf = new EncodingTreeNode {
        text[index], nullptr, nullptr
    };

    /* Last character? If so, that's all. */
    if (index + 1 == text.size()) return leaf;

    /* Otherwise, build a larger tree. */
    else return new EncodingTreeNode {
        ' ', leaf, strandTreeFor(text, index + 1)
    };
}

TEST("huffmanTreeFor reports errors on invalid inputs.") {
    EXPECT_ERROR(huffmanTreeFor(""));    // No characters
    EXPECT_ERROR(huffmanTreeFor("a"));   // One character
    EXPECT_ERROR(huffmanTreeFor("aaa")); // One character
}

TEST("huffmanTreeFor builds tree for two characters.") {
    EncodingTreeNode* reference = new EncodingTreeNode {
        ' ', new EncodingTreeNode {'a', nullptr, nullptr}, new EncodingTreeNode {'b', nullptr, nullptr}
    };

    EncodingTreeNode* tree = huffmanTreeFor("aaabbbb");
    EXPECT(isEncodingTree(tree));
    EXPECT(areEqual(tree, reference));

    deleteTree(reference);
    deleteTree(tree);
}

TEST("huffmanTreeFor works on the full range of characters.") {
    /* Get a string of all possible characters, then pair them off and see what we find. */
    string allChars = pangrammaticString();
    for (size_t i = 0; i < allChars.size(); i += 2) {
        string toEncode;
        toEncode += allChars[i];
        toEncode += allChars[i + 1];
        toEncode += allChars[i + 1];

        EncodingTreeNode* reference = new EncodingTreeNode {
            ' ',
            new EncodingTreeNode {allChars[i], nullptr, nullptr},
            new EncodingTreeNode {allChars[i + 1], nullptr, nullptr}
        };

        EncodingTreeNode* tree = huffmanTreeFor(toEncode);
        EXPECT(isEncodingTree(tree));
        EXPECT(areEqual(tree, reference));

        deleteTree(reference);
        deleteTree(tree);
    }
}



TEST("huffmanTreeFor uses cumulative weights (v1).") {
    /* This tree:
     *                 *
     *                / \
     *               *   D
     *              / \
     *             C   *
     *                / \
     *               A   B
     */
    EncodingTreeNode* reference = new EncodingTreeNode {
        '*',
            new EncodingTreeNode { '!',
                new EncodingTreeNode { 'C', nullptr, nullptr },
                new EncodingTreeNode { '?',
                    new EncodingTreeNode { 'A', nullptr, nullptr },
                    new EncodingTreeNode { 'B', nullptr, nullptr }
                }
            },
            new EncodingTreeNode { 'D', nullptr, nullptr }
    };

    /* Ax2, Bx3, Cx4, Dx10 */
    EncodingTreeNode* tree = huffmanTreeFor("AABBBCCCCDDDDDDDDDD");
    EXPECT(isEncodingTree(tree));
    EXPECT(areEqual(tree, reference));

    deleteTree(reference);
    deleteTree(tree);
}

TEST("huffmanTreeFor uses cumulative weights (v2).") {
    /*
     *          *
     *       /     \
     *      *       *
     *     / \     / \
     *    D   E   F   *
     *               / \
     *              C   *
     *                 / \
     *                A   B
     */
    EncodingTreeNode* reference =new EncodingTreeNode {
        ' ',
        new EncodingTreeNode {
            ' ',
            new EncodingTreeNode{ 'D', nullptr, nullptr },
            new EncodingTreeNode{ 'E', nullptr, nullptr }
        },
        new EncodingTreeNode {
            ' ',
            new EncodingTreeNode{ 'F', nullptr, nullptr },
            new EncodingTreeNode {
                ' ',
                new EncodingTreeNode{ 'C', nullptr, nullptr },
                new EncodingTreeNode{
                    ' ',
                    new EncodingTreeNode{ 'A', nullptr, nullptr },
                    new EncodingTreeNode{ 'B', nullptr, nullptr },
                }
            }
        }
    };

    /* Ax2, Bx3, Cx4, Dx6, Ex7, Fx8 */
    EncodingTreeNode* tree = huffmanTreeFor("AABBBCCCCDDDDDDEEEEEEEFFFFFFFF");
    EXPECT(isEncodingTree(tree));

    EXPECT(areEqual(tree, reference));

    deleteTree(reference);
    deleteTree(tree);
}

TEST("decodeText works on small sample.") {
    /* This tree:
     *                 *
     *                / \
     *               O   *
     *                  / \
     *                 *   N
     *                / \
     *               M   S
     */
    EncodingTreeNode* tree = new EncodingTreeNode {
        '*',
            new EncodingTreeNode { 'O', nullptr, nullptr },
            new EncodingTreeNode { '*',
                new EncodingTreeNode{ '*',
                    new EncodingTreeNode { 'M', nullptr, nullptr },
                    new EncodingTreeNode { 'S', nullptr, nullptr }
                },
                new EncodingTreeNode{ 'N', nullptr, nullptr }
            }
    };

    /* What you get if you encode MONSOON with this tree. */
    Queue<Bit> bits = { 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1 };

    EXPECT_EQUAL(decodeText(bits, tree), "MONSOON");

    deleteTree(tree);
}

TEST("Can decode all char values.") {
    /* All possible characters. */
    string allChars = pangrammaticString();

    /* Try decoding all pairs of adjacent characters. */
    for (size_t i = 0; i < allChars.size(); i += 2) {
        string expected;
        expected += allChars[i];
        expected += allChars[i + 1];
        expected += allChars[i + 1];

        EncodingTreeNode* tree = new EncodingTreeNode {
            ' ',
            new EncodingTreeNode {allChars[i], nullptr, nullptr},
            new EncodingTreeNode {allChars[i + 1], nullptr, nullptr}
        };

        /* Decode the bitstream 011, which should map back to the expected
         * string.
         */
        Queue<Bit> bits = { 0, 1, 1 };
        EXPECT_EQUAL(decodeText(bits, tree), expected);

        deleteTree(tree);
    }
}

PROVIDED_TEST("encodeText works on small sample.") {
    /* This tree:
     *                 *
     *                / \
     *               O   *
     *                  / \
     *                 *   N
     *                / \
     *               M   S
     */
    EncodingTreeNode* tree = new EncodingTreeNode {
           '*',
           new EncodingTreeNode { 'O', nullptr, nullptr },
               new EncodingTreeNode { '*',
               new EncodingTreeNode{ '*',
               new EncodingTreeNode { 'M', nullptr, nullptr },
               new EncodingTreeNode { 'S', nullptr, nullptr }
            },
            new EncodingTreeNode{ 'N', nullptr, nullptr }
        }
    };

    /* What you get if you encode MONSOON with this tree. */
    Queue<Bit> expected = { 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1 };

    EXPECT_EQUAL(encodeText("MONSOON", tree), expected);

    deleteTree(tree);
}

PROVIDED_TEST("Can encode all char values.") {
    /* All possible characters. */
    string allChars = pangrammaticString();

    /* Try encoding all pairs of adjacent characters. */
    for (size_t i = 0; i < allChars.size(); i += 2) {
        string toEncode;
        toEncode += allChars[i];
        toEncode += allChars[i + 1];
        toEncode += allChars[i + 1];

        EncodingTreeNode* tree = new EncodingTreeNode {
                ' ',
                new EncodingTreeNode {allChars[i], nullptr, nullptr},
                new EncodingTreeNode {allChars[i + 1], nullptr, nullptr}
        };

        /* See what bits we get back. We should get 011, since the first
         * character has code 0 and the second has code 1.
         */
        Queue<Bit> bits = encodeText(toEncode, tree);
        Queue<Bit> expected = { 0, 1, 1 };

        EXPECT_EQUAL(bits, expected);

        deleteTree(tree);
    }
}

PROVIDED_TEST("decodeText undoes encodeText on range of sample strings.") {
    Vector<string> testCases = {
        "THAT THAT IS IS THAT THAT IS NOT IS NOT IS THAT IT IT IS",
        "AABAAABBABAAABAAAA",
        ":-) :-D XD <(^_^)>",
        pangrammaticString(),
    };

    for (string test: testCases) {
        /* Use a silly encoding scheme, but one that works regardless of what
         * characters are provided.
         */
        EncodingTreeNode* tree = strandTreeFor(test);
        EXPECT(isEncodingTree(tree));

        Queue<Bit> bits = encodeText(test, tree);
        string result = decodeText(bits, tree);

        EXPECT_EQUAL(test.size(), result.size());
        EXPECT_EQUAL(test, result);

        deleteTree(tree);
    }
}

PROVIDED_TEST("Can decode an example tree.") {
    /* This encodes this tree:
     *
     *                 *
     *                / \
     *               *   C
     *              / \
     *             A   B
     */
    Queue<Bit>  bits   = { 1, 1, 0, 0, 0 };
    Queue<char> leaves = { 'A', 'B', 'C' };

    EncodingTreeNode* tree = decodeTree(bits, leaves);
    EXPECT(isEncodingTree(tree));

    /* Confirm this is the right tree. */
    EncodingTreeNode* expected = new EncodingTreeNode {
        '*',
            new EncodingTreeNode {
                '*',
                new EncodingTreeNode { 'A', nullptr, nullptr },
                new EncodingTreeNode { 'B', nullptr, nullptr },
            },
            new EncodingTreeNode { 'C', nullptr, nullptr }
    };

    EXPECT(areEqual(tree, expected));

    deleteTree(tree);
    deleteTree(expected);
}

PROVIDED_TEST("Can decode trees using all possible char values.") {
    /* All possible characters. */
    string allChars = pangrammaticString();

    /* Try encoding all pairs of adjacent characters. */
    for (size_t i = 0; i < allChars.size(); i += 2) {
        EncodingTreeNode* expected = new EncodingTreeNode {
            ' ',
            new EncodingTreeNode {allChars[i], nullptr, nullptr},
            new EncodingTreeNode {allChars[i + 1], nullptr, nullptr}
        };
        Queue<Bit>  treeBits   = { 1, 0, 0 };
        Queue<char> treeLeaves = { allChars[i], allChars[i + 1] };

        EncodingTreeNode* tree = decodeTree(treeBits, treeLeaves);
        EXPECT(isEncodingTree(tree));
        EXPECT(areEqual(tree, expected));

        deleteTree(tree);
        deleteTree(expected);
    }
}

PROVIDED_TEST("Can encode an example tree.") {
    /* Build an encoding tree for "ABBCCCC." It should look like this:
     *
     *                 *
     *                / \
     *               *   C
     *              / \
     *             A   B
     *
     * This will compress down to
     *
     *        11000
     *        ABC
     */
    EncodingTreeNode* tree = huffmanTreeFor("ABBCCCC");

    Queue<Bit>  bits;
    Queue<char> leaves;

    encodeTree(tree, bits, leaves);

    Queue<Bit>  expectedBits   = { 1, 1, 0, 0, 0 };
    Queue<char> expectedLeaves = { 'A', 'B', 'C' };

    EXPECT_EQUAL(bits,   expectedBits);
    EXPECT_EQUAL(leaves, expectedLeaves);

    deleteTree(tree);
}

PROVIDED_TEST("Can encode trees using all possible char values.") {
    /* All possible characters. */
    string allChars = pangrammaticString();

    /* Try encoding all pairs of adjacent characters. */
    for (size_t i = 0; i < allChars.size(); i += 2) {
        EncodingTreeNode* tree = new EncodingTreeNode {
            ' ',
            new EncodingTreeNode {allChars[i], nullptr, nullptr},
            new EncodingTreeNode {allChars[i + 1], nullptr, nullptr}
        };

        /* See what we get back. It should be the bitstring 100 (root with
         * two children) and the two leaves, in order.
         */
        Queue<Bit>  treeBits;
        Queue<char> treeLeaves;

        Queue<Bit>  expectedBits = { 1, 0, 0 };
        Queue<char> expectedLeaves = { allChars[i], allChars[i + 1] };

        encodeTree(tree, treeBits, treeLeaves);
        EXPECT_EQUAL(treeBits, expectedBits);
        EXPECT_EQUAL(treeLeaves, expectedLeaves);

        deleteTree(tree);
    }
}

PROVIDED_TEST("decodeTree undoes encodeTree on sample strings.") {
    /* Make a Huffman tree for the string of all characters. */
    EncodingTreeNode* sourceTree = huffmanTreeFor(pangrammaticString());
    EXPECT(isEncodingTree(sourceTree));

    /* Encode, then decode it. */
    Queue<Bit>  bits;
    Queue<char> leaves;
    encodeTree(sourceTree, bits, leaves);

    EncodingTreeNode* resultTree = decodeTree(bits, leaves);
    EXPECT(isEncodingTree(resultTree));
    EXPECT(areEqual(sourceTree, resultTree));

    deleteTree(sourceTree);
    deleteTree(resultTree);
}

PROVIDED_TEST("Can decompress a small sample file.") {
    HuffmanResult file = {
        { 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0 },
        { 'u', 'k', 'p', 'n', 'a', 'm', 'h' },
        { 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1,
          0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0,
          0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0 }
    };

    EXPECT_EQUAL(decompress(file), "humuhumunukunukuapuaa");
}

PROVIDED_TEST("Compress reports errors on bad inputs.") {
    EXPECT_ERROR(compress(""));
    EXPECT_ERROR(compress("A"));
    EXPECT_ERROR(compress("AAAA"));
}

PROVIDED_TEST("Can compress a small sample file.") {
    HuffmanResult file = compress("ABANANAABANDANA");
    Queue<Bit>  treeBits    = { 1, 1, 1, 0, 0, 0, 0 };
    Queue<char> treeChars   = { 'D', 'B', 'N', 'A' };
    Queue<Bit>  messageBits = { 1, 0, 0, 1, 1, 0, 1, 1, 0,
                                1, 1, 1, 0, 0, 1, 1, 0, 1,
                                0, 0, 0, 1, 0, 1, 1 };

    EXPECT_EQUAL(file.treeBits, treeBits);
    EXPECT_EQUAL(file.treeLeaves, treeChars);
    EXPECT_EQUAL(file.messageBits, messageBits);
}

PROVIDED_TEST("Compress undoes decompress on a range of strings.") {
    Vector<string> testCases = {
        "THAT THAT IS IS THAT THAT IS NOT IS NOT IS THAT IT IT IS",
        "AABAAABBABAAABAAAA",
        ":-) :-D XD <(^_^)>",
        pangrammaticString(),
    };

    for (string test: testCases) {
        HuffmanResult file = compress(test);
        string result = decompress(file);

        EXPECT_EQUAL(result.size(), test.size());
        EXPECT_EQUAL(test, result);
    }
}
