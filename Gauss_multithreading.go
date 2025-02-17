package main

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"strconv"
	"strings"
	"sync" // Пакет sync для синхронизации горутин
	"time"
)

// Функция решает систему линейных алгебраических уравнений (СЛАУ) методом Гаусса
func solveSLAUParallelGauss(matrix [][]float64) ([]float64, error) {
	n := len(matrix)        // количество строк (уравнений)
	m := len(matrix[0]) - 1 // количество столбцов коэффициентов (переменных)

	if n != m {
		return nil, fmt.Errorf("количество уравнений не равно количеству переменных, система не квадратная")
	}

	startTime := time.Now() // Запускаем таймер

	// Прямой ход метода Гаусса
	for i := 0; i < n; i++ {
		pivotRow := i
		for j := i + 1; j < n; j++ {
			if math.Abs(matrix[j][i]) > math.Abs(matrix[pivotRow][i]) {
				pivotRow = j
			}
		}

		if pivotRow != i {
			matrix[i], matrix[pivotRow] = matrix[pivotRow], matrix[i]
		}

		if math.Abs(matrix[i][i]) < 1e-9 {
			allZerosRow := true
			for k := 0; k < m; k++ {
				if math.Abs(matrix[i][k]) > 1e-9 {
					allZerosRow = false
					break
				}
			}
			if allZerosRow && math.Abs(matrix[i][m]) > 1e-9 {
				return nil, fmt.Errorf("система несовместна, нет решений")
			}
			if math.Abs(matrix[i][i]) < 1e-9 {
				continue
			}
		}

		var wg sync.WaitGroup // Объявляем WaitGroup для ожидания завершения всех горутин
		// WaitGroup используется для синхронизации(основная горутина будет ждать,
		// пока все запущенные на параллельное выполнение горутины не завершат свою работу.)

		for j := i + 1; j < n; j++ {
			wg.Add(1)                 // Увеличиваем счетчик WaitGroup на 1 перед запуском каждой горутины
			go func(currentRow int) { // Запускаем горутину
				defer wg.Done() // Уменьшаем счетчик WaitGroup на 1 после завершения работы горутины
				factor := matrix[currentRow][i] / matrix[i][i]
				for k := i; k <= m; k++ {
					matrix[currentRow][k] -= factor * matrix[i][k]
				}
			}(j)
		}
		wg.Wait()
	}

	// Обратный ход метода Гаусса
	solution := make([]float64, n)
	for i := n - 1; i >= 0; i-- {
		sum := 0.0
		for j := i + 1; j < n; j++ {
			sum += matrix[i][j] * solution[j]
		}
		solution[i] = (matrix[i][m] - sum) / matrix[i][i]
	}

	elapsedTime := time.Since(startTime)
	fmt.Printf("Время вычисления: %s\n", elapsedTime)

	return solution, nil
}

// Функция читает расширенную матрицу из файла.
// Файл должен содержать строки, представляющие строки матрицы, числа в строке должны быть разделены пробелами
func readMatrixFromFile(filename string) ([][]float64, error) {
	file, err := os.Open(filename) // Открываем файл для чтения
	if err != nil {
		return nil, fmt.Errorf("не удалось открыть файл: %w", err) // Возвращаем ошибку, если не удалось открыть файл
	}
	defer file.Close() // Убеждаемся, что файл будет закрыт после завершения функции

	var matrix [][]float64            // Инициализируем матрицу
	scanner := bufio.NewScanner(file) // Создаем сканер для чтения файла построчно

	rowNum := 0          // Счетчик строк для сообщений об ошибках
	for scanner.Scan() { // Читаем файл построчно
		rowNum++
		line := scanner.Text()         // Получаем текущую строку
		fields := strings.Fields(line) // Разбиваем строку на поля (числа, разделенные пробелами)

		var row []float64              // Инициализируем строку матрицы
		for _, field := range fields { // Итерируем по полям в строке
			value, err := strconv.ParseFloat(field, 64) // Пытаемся преобразовать поле в число float64
			if err != nil {
				return nil, fmt.Errorf("ошибка парсинга числа '%s' в строке %d: %w", field, rowNum, err) // Возвращаем ошибку, если не удалось преобразовать
			}
			row = append(row, value) // Добавляем число в текущую строку
		}
		matrix = append(matrix, row) // Добавляем строку в матрицу
	}

	if err := scanner.Err(); err != nil { // Проверяем наличие ошибок сканирования файла
		return nil, fmt.Errorf("ошибка чтения файла: %w", err) // Возвращаем ошибку, если произошла ошибка чтения
	}

	if len(matrix) == 0 {
		return nil, fmt.Errorf("файл matrix.txt пуст или не содержит данных матрицы")
	}

	// Проверка на то, что все строки матрицы имеют одинаковое количество столбцов
	cols := len(matrix[0])
	for i := 1; i < len(matrix); i++ {
		if len(matrix[i]) != cols {
			return nil, fmt.Errorf("строка %d имеет другое количество столбцов, чем первая строка", i+1)
		}
	}

	return matrix, nil // Возвращаем матрицу и nil, если успешно прочитано
}

func main() {
	filename := "matrix.txt"                    // Имя файла с матрицей
	matrix, err := readMatrixFromFile(filename) // Читаем матрицу из файла
	if err != nil {
		fmt.Println("Ошибка при чтении матрицы из файла:", err) // Выводим сообщение об ошибке, если есть
		return                                                  // Завершаем программу с ошибкой
	}

	solution, err := solveSLAUParallelGauss(matrix) // Решаем СЛАУ, используя матрицу из файла
	if err != nil {
		fmt.Println("Ошибка при решении СЛАУ:", err) // Выводим сообщение об ошибке, если есть
		return                                       // Завершаем программу с ошибкой
	}

	fmt.Println("Решение СЛАУ:") // Выводим результаты
	for i, sol := range solution {
		fmt.Printf("x%d = %f\n", i+1, sol) // Выводим каждое решение
	}
}
