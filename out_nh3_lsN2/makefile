CC = gcc
CFLAGS = -O2 -Wall -Wextra -I. -Ijacobs
LDFLAGS = -llapack -lblas -lm  # 数学ライブラリをリンク

# メインディレクトリのソースファイル
SRC = main.c chem_utils.c dydt.c jacob.c mass_mole.c mechanism.c rxn_rates.c rxn_rates_pres_mod.c sparse_multiplier.c spec_rates.c
OBJ = $(SRC:.c=.o)

# jacobs ディレクトリ内のソースファイルを取得
JACOBS_SRC = $(wildcard jacobs/*.c)
JACOBS_OBJ = $(JACOBS_SRC:.c=.o)

# すべてのオブジェクトファイル
OBJ += $(JACOBS_OBJ)

# 依存関係
DEPS = chem_utils.h dydt.h jacob.h mass_mole.h mechanism.h rates.h rxn_rates.h sparse_multiplier.h spec_rates.h header.h $(wildcard jacobs/*.h)

# 実行ファイルの名前
TARGET = main

# すべてのオブジェクトファイルをリンクして実行ファイルを作成
$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# 各Cファイルをコンパイル
%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<

# クリーンアップ
.PHONY: clean
clean:
	rm -f $(OBJ) $(TARGET)
