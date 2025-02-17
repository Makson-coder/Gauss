package main

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
	"strings"
)

// Функция для чтения расширенной матрицы из файла
func readMatrixFromFile(filename string) ([][]float64, error) {
	file, err := os.Open(filename)
	if err != nil { // Проверяем, произошла ли ошибка при открытии файла.
		return nil, fmt.Errorf("ошибка открытия файла: %w", err)
	}
	defer file.Close()

	var matrix [][]float64 //Создаем двумерный слайс float64 для хранения матрицы
	scanner := bufio.NewScanner(file)
	for scanner.Scan() { //Запускаем цикл для чтения файла построчно.
		line := scanner.Text()
		if line == "" { // Пропускаем пустые строки
			continue //Если строка пустая, переходим к следующей итерации цикла
		}
		strValues := strings.Fields(line) // Разделяем текущую строку "line" на подстроки (слова), используя пробелы в качестве разделителей
		row := make([]float64, len(strValues))
		for i, strVal := range strValues {
			val, err := strconv.ParseFloat(strVal, 64) //Пытаемся преобразовать текущую подстроку "strVal" в число типа float64
			if err != nil {
				return nil, fmt.Errorf("ошибка парсинга числа '%s': %w", strVal, err)
			}
			row[i] = val
		}
		matrix = append(matrix, row) // Добавляем сформированную строку чисел "row" в общую матрицу "matrix"
	}

	if err := scanner.Err(); err != nil { // Проверяем, не возникло ли ошибок при сканировании файла после завершения цикла чтения строк.
		return nil, fmt.Errorf("ошибка чтения файла: %w", err)
	}

	if len(matrix) == 0 { // Проверяем, является ли матрица пустой после чтения из файла.
		return nil, fmt.Errorf("файл матрицы пуст или не содержит корректных данных")
	}
	if len(matrix) != len(matrix[0])-1 {
		return nil, fmt.Errorf("матрица должна быть квадратной для метода Гаусса и иметь столбец свободных членов. Количество строк должно быть равно количеству столбцов коэффициентов.")
	}

	return matrix, nil
}

// Функция для решения СЛАУ методом Гаусса и проверки совместности
func gaussianElimination(matrix [][]float64) ([]float64, error) {
	n := len(matrix) // Получаем количество строк в матрице, которое соответствует размерности системы уравнений (n x n).
	if n == 0 {
		return nil, fmt.Errorf("матрица не может быть пустой")
	}

	a := make([][]float64, n) // Копия матрицы, чтобы не изменять оригинал
	for i := 0; i < n; i++ {
		a[i] = make([]float64, len(matrix[i]))
		copy(a[i], matrix[i])
	}

	for i := 0; i < n; i++ {
		// Поиск главного элемента в i-м столбце
		pivot := i 
		for j := i + 1; j < n; j++ {
			if abs(a[j][i]) > abs(a[pivot][i]) {
				pivot = j
			}
		}

		// Перестановка строк, чтобы главный элемент был на диагонали
		a[i], a[pivot] = a[pivot], a[i]

		if a[i][i] == 0 {
			// Проверка на несовместность: если главный элемент 0, проверяем всю строку
			allZerosRow := true
			for k := i; k < n; k++ { // Проверяем только коэффициенты
				if abs(a[i][k]) > 1e-9 {
					allZerosRow = false
					break
				}
			}
			if allZerosRow && abs(a[i][n]) > 1e-9 { // Если все коэффициенты 0, а свободный член нет
				return nil, fmt.Errorf("система несовместна, решения не существует")
			}
			// Если главный элемент 0, но строка не ведет к несовместности, продолжаем (может быть бесконечно много решений или вырожденная матрица)
			// В данном коде для квадратной системы вырожденная матрица тоже приведет к ошибке, но другого рода (деление на 0 на следующих шагах, если не обработать специально)
			if i == n-1 && abs(a[i][i]) < 1e-9 && abs(a[i][n]) > 1e-9 {
				return nil, fmt.Errorf("система несовместна, решения не существует (последняя строка)")
			}
			if i < n-1 { // Если не последняя строка, просто пропускаем нулевой главный элемент и продолжаем, это может быть случай бесконечного числа решений, но для квадратной системы - скорее всего вырожденность
				continue // для простоты в данном коде, если встретили нулевой главный элемент не в последней строке, просто пропускаем, что может привести к делению на ноль далее или неверному результату для вырожденных систем
			}

		}

		// Прямой ход (исключение Гаусса)
		for j := i + 1; j < n; j++ {
			factor := a[j][i] / a[i][i]
			for k := i; k <= n; k++ {
				a[j][k] -= factor * a[i][k]
			}
		}
	}

	// Проверка на несовместность после прямого хода (на случай, если не обнаружили раньше)
	for i := 0; i < n; i++ {
		allZerosRow := true
		for j := 0; j < n; j++ {
			if abs(a[i][j]) > 1e-9 {
				allZerosRow = false
				break
			}
		}
		if allZerosRow && abs(a[i][n]) > 1e-9 {
			return nil, fmt.Errorf("система несовместна, решения не существует (после прямого хода)")
		}
	}

	// Обратный ход (обратная подстановка)
	solution := make([]float64, n)
	for i := n - 1; i >= 0; i-- {
		if abs(a[i][i]) < 1e-9 { // Обработка случая нулевого главного элемента после преобразований (для случая вырожденности или бесконечного числа решений)
			if abs(a[i][n]) > 1e-9 { // Если на диагонали 0, а в столбце свободных членов не 0, то система несовместна (хотя проверка выше должна была это поймать)
				return nil, fmt.Errorf("система несовместна, решения не существует (обратный ход, нулевой главный элемент)")
			}
			solution[i] = 0 // Можно присвоить 0 или NaN, если хотим явно указать на неопределенность, здесь для простоты 0,  но для вырожденной системы может быть бесконечно много решений, что данный код не обрабатывает корректно
		} else {
			sum := 0.0
			for j := i + 1; j < n; j++ {
				sum += a[i][j] * solution[j]
			}
			solution[i] = (a[i][n] - sum) / a[i][i]
		}
	}

	return solution, nil
}

// Функция для проверки решения СЛАУ
func checkSolution(matrix [][]float64, solution []float64) bool {
	n := len(matrix)
	for i := 0; i < n; i++ {
		sum := 0.0
		for j := 0; j < n; j++ {
			sum += matrix[i][j] * solution[j]
		}
		if abs(sum-matrix[i][n]) > 1e-9 { // Используем небольшую погрешность для сравнения float
			fmt.Printf("Уравнение %d не выполняется. Левая часть: %f, Правая часть: %f\n", i+1, sum, matrix[i][n]) // Выводим информацию о проблеме
			return false
		}
	}
	return true
}

// Вспомогательная функция для вычисления абсолютного значения float64
func abs(x float64) float64 {
	if x < 0 {
		return -x
	}
	return x
}

func main() {
	filename := "matrix.txt" // Жестко заданное имя файла матрицы
	matrix, err := readMatrixFromFile(filename)
	if err != nil {
		fmt.Println("Ошибка при чтении матрицы из файла:", err)
		return
	}

	fmt.Println("Исходная расширенная матрица:")
	for _, row := range matrix {
		fmt.Println(row)
	}

	solution, err := gaussianElimination(matrix)
	if err != nil {
		fmt.Println("Ошибка при решении СЛАУ методом Гаусса:", err)
		if strings.Contains(err.Error(), "несовместна") {
			fmt.Println("Система несовместна, решения не существует.")
		}
		return
	}

	fmt.Println("\nРешение СЛАУ:")
	for i, val := range solution {
		fmt.Printf("x%d = %f\n", i+1, val)
	}

	if checkSolution(matrix, solution) {
		fmt.Println("\nРешение СЛАУ корректно.")
	} else {
		fmt.Println("\nРешение СЛАУ НЕ корректно или не является точным из-за вычислительных погрешностей.")
	}
}
